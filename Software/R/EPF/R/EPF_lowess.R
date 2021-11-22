EPF_lowess <- function(data_ddf, var_x, var_y, x_0, var_x_range = NULL, alpha = 0.01, K = 2, J = c(4, 8, 8), include_residuals = TRUE){
  
  data2var_x_scale <- function(data_ddf, var_x){
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_x]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      collect("minimum", -min(data_subset))
      collect("maximum", max(data_subset))
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- max(output, unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_x = var_x)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys <- unlist(lapply(mr_out, function(l) l[[1]]))
    values <- unlist(lapply(mr_out, "[[", 2))
    
    minimum <- -values[keys == "minimum"]
    maximum <- values[keys == "maximum"]
    
    M = 1/(maximum-minimum)
    C = -minimum/(maximum-minimum)
    
    var_x_scale = list(M = M, C = C)        
    
    return(var_x_scale)
  }
  
  data2lowess_WEP_3_parameter <- function(data_ddf, var_x, J_1, J_2, J_3, M, C){
    K = 2*J_1*J_2*J_3
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_x]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      scaled_data_subset <- M*as.vector(data_subset) + C
      T <- get_trigonometric_summand_3_parameter(scaled_data_subset, J_1, J_2, J_3)
      for (i in 1:length(T)){collect(i, T[i])}
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- output + sum(unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_x = var_x, M = M, C = C, J_1 = J_1, J_2 = J_2, J_3 = J_3)
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
  
  data2lowess_WEP_4_parameter <- function(data_ddf, var_x, J_1, J_2, J_3, J_4, M, C){
    K = 2*J_1*J_2*J_3*J_4
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_x]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      scaled_data_subset <- M*as.vector(data_subset) + C
      T <- get_trigonometric_summand_4_parameter(scaled_data_subset, J_1, J_2, J_3, J_4)
      for (i in 1:length(T)){collect(i, T[i])}
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- output + sum(unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_x = var_x, M = M, C = C, J_1 = J_1, J_2 = J_2, J_3 = J_3, J_4 = J_4)
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
  
  Trig_WEP2lowess_widths <- function(Trig_WEP, x_0, alpha, M){
    J <- as.integer((length(Trig_WEP)-1)/2)
    h <- rep(0, length(x_0))
    x_0_scaled <- M*x_0+C
    
    for (i in 1:length(x_0)){
      h[i]=uniroot(fn <- function(m){return(sum(Trig_WEP*get_lowess_multiplier(x_0_scaled[i], m, J))-alpha)}, c(0, 1))$root
    }
    
    lowess_widths <- h/M
    return(lowess_widths)
  }
  
  data2lowess_predicted_values <- function(data_ddf, var_x, var_y, x_0, h, K){
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, c(var_x, var_y)]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      P <- get_loess_powers(data_subset[, var_x], data_subset[, var_y], x_0, h, K)
      for (i in 1:nrow(P)){
        for (k in 1:ncol(P)){
          collect(c(i, k), P[i, k])
        }
      }
    }
    })
    reduce <- expression(pre = {sum <- 0}, reduce = {sum <- sum + sum(unlist(reduce.values))}, post = {collect(reduce.key, sum)})
    parameter_list <- list(var_x = var_x, var_y = var_y, x_0 = x_0, h = h, K = K)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    keys_id <- sapply(lapply(mr_out, function(x) x[[1]][1]), as.integer)
    keys_power <- sapply(lapply(mr_out, function(x) x[[1]][2]), as.integer)
    values <- unlist(lapply(mr_out, "[[", 2))
    
    lowess_power_matrix <- matrix(0, length(x_0), 3*K+2)
    for (i in 1:length(x_0)){
      for (k in 1:(3*K+2)){
        lowess_power_matrix[i, k] = values[keys_id == i & keys_power == k]
      }
    }
    
    y_0_hat <- apply(lowess_power_matrix, 1, lowess_powers2predicted_value)
    
    return(y_0_hat)
  }

  lowess_powers2predicted_value <- function(p){
    K <- (length(p)-2)/3
    xx_power <- matrix(0, K+1, K+1)
    
    for(k in 1:(K+1)){
      for(l in 1:(K+1)){
        xx_power[k, l]=p[k+l-1]
      }
    }
    
    xy_power <- p[(2*K+2):(3*K+2)]
    y_0_vector <- solve(xx_power)%*%xy_power
    y_0_hat <- y_0_vector[1]
    
    return(y_0_hat)
  }
  
  data2residuals <- function(data_ddf, var_x, var_y, x_0, y_0_hat){
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, c(var_x, var_y)]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (nrow(data_subset) > 0){
      for(i in 1:length(x_0)){
        T = c()
        for(n in 1:nrow(data_subset)){
          if(data_subset[n, var_x] == x_0[i]){
            T = c(T, data_subset[n, var_y])
          }
        }
        collect(i, T)
      }
    }
    })
    reduce <- expression(pre = {output <- c()}, reduce = {output <- c(output, unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_x = var_x, var_y = var_y, x_0 = x_0)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys <- unlist(lapply(mr_out, function(x) x[[1]]))
    values <- lapply(mr_out, "[[", 2)

    y_actual = list()

    for(i in 1:length(x_0)){
      y_actual[[keys[i]]] = values[[i]]
    }

    lowess_stats = list()

    for(i in 1:length(x_0)){
      if(!is.null(values[[i]])){
        lowess_stats[[i]] = list(x = x_0[i], y_actual_values = y_actual[[i]], residuals = y_actual[[i]]-y_0_hat[i])
      }
      else{
        lowess_stats[[i]] = list(x = x_0[i], y_actual_values = NULL, residuals = NULL)
      }
    }
    
    return(lowess_stats)
  }
    
  start_time <- proc.time()
  
  if(!inherits(data_ddf, "ddf")) stop("data_ddf should be a distributed data frame.")
  
  if(!is.character(var_x) | var_x %in% names(attributes(data_ddf)$ddf$vars) == FALSE) stop("var_x should be a character representing the expalanatory variable present in the data.")
  
  if(!is.character(var_y) | var_y %in% names(attributes(data_ddf)$ddf$vars) == FALSE) stop("var_y should be a character representing the response variable present in the data.")
  
  if(is.null(var_x_range) == TRUE){
    var_x_scale = data2var_x_scale(data_ddf, var_x)
    M = var_x_scale$M
    C = var_x_scale$C
    var_x_range = c((-C)/M, (1-C)/M)
  }
  else{
    M = 1/(var_x_range[2]-var_x_range[1])
    C = -var_x_range[1]/(var_x_range[2]-var_range[1])
  }
  
  if(!is.vector(x_0) | !is.numeric(x_0) | min(x_0) < var_x_range[1] | max(x_0) > var_x_range[2]) stop("x_0 should be a vector of test-data and it should have the same range as var_x")
  
  if(!is.double(alpha) | alpha <= 0 | alpha >= 1) stop("alpha should be fraction in (0, 1)")
  
  if(!is.double(K) | !K == floor(K)) stop("K should be an integer")
  
  if(!is.vector(J) | !is.numeric(J) | !all(J==floor(J)) | length(J) %in% c(3, 4) == FALSE) stop("J should be and integer vector containing 3 or 4 parameters.")
  
  if(length(J) == 3) Trig_WEP = data2lowess_WEP_3_parameter(data_ddf, var_x, J_1 = J[1], J_2 = J[2], J_3 = J[3], M, C)
  
  if(length(J) == 4) Trig_WEP = data2lowess_WEP_4_parameter(data_ddf, var_x, J_1 = J[1], J_2 = J[2], J_3 = J[3], J_4 = J[4], M, C)
  
  h = Trig_WEP2lowess_widths(Trig_WEP, x_0, alpha, M)
  y_0_hat = data2lowess_predicted_values(data_ddf, var_x, var_y, x_0, h, K)
  
  if(include_residuals == TRUE){
    lowess_statistics = data2residuals(data_ddf, var_x, var_y, x_0, y_0_hat)
    error = unlist(lapply(lowess_statistics, function(x) x[[3]]))
    RSS = sum(error^2)
  }
  
  stop_time = proc.time()
  
  if(include_residuals == TRUE){
    EPF_lowess_output = list(Lowess_fit = data.frame(x = x_0, y_hat = y_0_hat), Lowess_statistics = lowess_statistics, Residual_sum_of_squares = RSS, Computation_time_in_seconds = unname(stop_time-start_time))
  } else {
    EPF_lowess_output = list(Lowess_fit = data.frame(x = x_0, y_hat = y_0_hat), Computation_time_in_seconds = unname(stop_time-start_time))
  }
  
  return(EPF_lowess_output)
}