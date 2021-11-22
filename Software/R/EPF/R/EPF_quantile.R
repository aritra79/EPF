EPF_quantile <- function(data_ddf, var_name, probs = seq(0, 1, 0.01), J = c(4, 8, 8), var_range = NULL){
  
  data2var_scale <- function(data_ddf, var_name){
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_name]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      collect("minimum", -min(data_subset))
      collect("maximum", max(data_subset))
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- max(output, unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_name = var_name)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys <- unlist(lapply(mr_out, function(l) l[[1]]))
    values <- unlist(lapply(mr_out, "[[", 2))
    
    minimum <- -values[keys == "minimum"]
    maximum <- values[keys == "maximum"]
    
    M = 1/(maximum-minimum)
    C = -minimum/(maximum-minimum)
    
    var_scale = list(M = M, C = C)        
    
    return(var_scale)
  }
  
  data2quantile_WEP_3_parameter <- function(data_ddf, var_name, J_1, J_2, J_3, M, C){
    K = 2*J_1*J_2*J_3
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_name]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      scaled_data_subset <- M*as.vector(data_subset) + C
      T <- get_trigonometric_summand_3_parameter(scaled_data_subset, J_1, J_2, J_3)
      for (i in 1:length(T)){collect(i, T[i])}
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- output + sum(unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_name = var_name, M = M, C = C, J_1 = J_1, J_2 = J_2, J_3 = J_3)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys <- unlist(lapply(mr_out, function(x) x[[1]][1]))
    values <- unlist(lapply(mr_out, "[[", 2))
    
    Trig_3D_power_SEP <- values[invPerm(as.integer(keys))]
    Trig_3D_power_WEP <- Trig_3D_power_SEP/Trig_3D_power_SEP[1]
    
    Trig_WEP <- Trig_3D_power_WEP2Trig_WEP(Trig_3D_power_WEP, J_1, J_2, J_3)
    
    return(Trig_WEP)
  }
  
  data2quantile_WEP_4_parameter <- function(data_ddf, var_name, J_1, J_2, J_3, J_4, M, C){
    K = 2*J_1*J_2*J_3*J_4
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_name]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      scaled_data_subset <- M*as.vector(data_subset) + C
      T <- get_trigonometric_summand_4_parameter(scaled_data_subset, J_1, J_2, J_3, J_4)
      for (i in 1:length(T)){collect(i, T[i])}
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- output + sum(unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_name = var_name, M = M, C = C, J_1 = J_1, J_2 = J_2, J_3 = J_3, J_4 = J_4)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys <- unlist(lapply(mr_out, function(x) x[[1]][1]))
    values <- unlist(lapply(mr_out, "[[", 2))
    
    Trig_4D_power_SEP <- values[invPerm(as.integer(keys))]
    Trig_4D_power_WEP <- Trig_4D_power_SEP/Trig_4D_power_SEP[1]
    
    Trig_WEP <- Trig_4D_power_WEP2Trig_WEP(Trig_4D_power_WEP, J_1, J_2, J_3, J_4)
    
    return(Trig_WEP)
  }
  
  Trig_3D_power_WEP2Trig_WEP <- function(Trig_3D_power_WEP, J_1, J_2, J_3) 
  {
    Trig_WEP_1 <- Trig_3D_power_WEP[1]
    Trig_3D_power_WEP <- matrix(Trig_3D_power_WEP[-1], nrow = 2)
    
    Cosine_3D_power_WEP <- array(Trig_3D_power_WEP[1, ], dim = c(J_1, J_2, J_3))
    Sine_3D_power_WEP <- array(Trig_3D_power_WEP[2, ], dim = c(J_1, J_2, J_3))
    
    Cosine_3D_product_WEP <- apply(apply(apply(Cosine_3D_power_WEP, c(1, 2), ChevTrB), c(1, 2), ChevTrB), c(1, 2), ChevTrA)
    Sine_3D_product_WEP <- apply(apply(apply(Sine_3D_power_WEP, c(1, 2), ChevTrB), c(1, 2), ChevTrB), c(1, 2), FUN = function(x) OddNegTr(ChevTrA(SineTr(x))))
    
    Cosine_WEP <- CosineSumTr(apply(Cosine_3D_product_WEP, 3, CosineSumTr))
    Sine_WEP <- SineSumTr(apply(Sine_3D_product_WEP, 3, SineSumTr))
    
    Trig_WEP <- c(Trig_WEP_1, rbind(Cosine_WEP, Sine_WEP))
    
    return(Trig_WEP)
  }
  
  Trig_4D_power_WEP2Trig_WEP <- function(Trig_4D_power_WEP, J_1, J_2, J_3, J_4) 
  {
    Trig_WEP_1 <- Trig_4D_power_WEP[1]
    Trig_4D_power_WEP <- matrix(Trig_4D_power_WEP[-1], nrow = 2)
    
    Cosine_4D_power_WEP <- array(Trig_4D_power_WEP[1, ], dim = c(J_1, J_2, J_3, J_4))
    Sine_4D_power_WEP <- array(Trig_4D_power_WEP[2, ], dim = c(J_1, J_2, J_3, J_4))
    
    Cosine_4D_product_WEP <- apply(apply(apply(apply(Cosine_4D_power_WEP, c(1, 2, 3), ChevTrB), c(1, 2, 3), ChevTrB), c(1, 2, 3), ChevTrB), c(1, 2, 3), ChevTrA)
    Sine_4D_product_WEP <- apply(apply(apply(apply(Sine_4D_power_WEP, c(1, 2, 3), ChevTrB), c(1, 2, 3), ChevTrB), c(1, 2, 3), ChevTrB), c(1, 2, 3), FUN = function(x) OddNegTr(ChevTrA(SineTr(x))))
    
    Cosine_WEP <- CosineSumTr(apply(apply(Cosine_4D_product_WEP, c(3, 4), CosineSumTr), 3, CosineSumTr))
    Sine_WEP <- SineSumTr(apply(apply(Sine_4D_product_WEP, c(3, 4), SineSumTr), 3, SineSumTr))
    
    Trig_WEP <- c(Trig_WEP_1, rbind(Cosine_WEP, Sine_WEP))
    
    return(Trig_WEP)
  }
  
  Trig_WEP2quantiles <- function(Trig_WEP, probs, M, C){
    J <- as.integer((length(Trig_WEP)-1)/2)
    Q <- rep(0, length(probs))

    for (i in 1:length(probs)) {
      if (probs[i]==0){
        Q[i] = -C/M
      } else if(probs[i]==1) {
        Q[i] = (1-C)/M
      } else {
        Q[i]=(optimize(fn <- function(m){return(sum(Trig_WEP*get_quantile_multiplier(m, J))+m*(probs[i]-0.5))}, interval = c(0, 1), maximum = TRUE)$maximum-C)/M
      }
    }

    return(Q)
  }
  
  start_time <- proc.time()
  
  if(!inherits(data_ddf, "ddf")) stop("data_ddf should be a distributed data frame.")
  
  if(!is.character(var_name) | var_name %in% names(attributes(data_ddf)$ddf$vars) == FALSE) stop("var should be a character representing a variable present in the data.")
  
  if(!is.vector(probs) | min(probs) < 0 | max(probs) > 1) stop("probs should be a vector of f-values")
  
  if(is.null(var_range) == TRUE){
    var_scale = data2var_scale(data_ddf, var_name)
    M = var_scale$M
    C = var_scale$C
    }
  else{
    M = 1/(var_range[2]-var_range[1])
    C = -var_range[1]/(var_range[2]-var_range[1])
  }
  
  if(!is.vector(J) | !all(J==floor(J)) | length(J) %in% c(3, 4) == FALSE) stop("J should be and integer vector containing 3 or 4 parameters.")
  
  if(length(J) == 3) Trig_WEP = data2quantile_WEP_3_parameter(data_ddf, var_name, J_1 = J[1], J_2 = J[2], J_3 = J[3], M, C)
  
  if(length(J) == 4) Trig_WEP = data2quantile_WEP_4_parameter(data_ddf, var_name, J_1 = J[1], J_2 = J[2], J_3 = J[3], J_4 = J[4], M, C)
  
  Q <- Trig_WEP2quantiles(Trig_WEP, probs, M, C)
  
  stop_time <- proc.time()
  
  EPF_quantile_output = list(Quantiles = data.frame(f_values = probs, Quantile_values = Q), Computation_time_in_seconds = unname(stop_time[3] - start_time[3]))
  
  return(EPF_quantile_output)
}
