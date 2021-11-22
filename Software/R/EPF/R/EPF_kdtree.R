EPF_kdtree <- function (data_ddf, var_names, D = 6, J = c(3, 4), var_ranges = NULL, include_cell_counts = TRUE){
  
  data2var_scales <- function(data_ddf, var_names){
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_names]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      min_subset <- unlist(lapply(data_subset, min))
      max_subset <- unlist(lapply(data_subset, max))
      for (v in c(var_names)){
        collect(c(v, "min"), -min_subset[v])
        collect(c(v, "max"), max_subset[v])
      }
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- max(output, unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_names = var_names)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys_var_names <- unlist(lapply(mr_out, function(l) l[[1]][1]))
    keys_type <- unlist(lapply(mr_out, function(l) l[[1]][2]))
    values <- unlist(lapply(mr_out, "[[", 2))
    
    var_mins <- c()
    var_maxs <- c()
    
    for (v in c(var_names)){
      var_mins <- c(var_mins, -values[keys_var_names == v & keys_type == "min"])
      var_maxs <- c(var_maxs, values[keys_var_names == v & keys_type == "max"])
    }
    
    M <- 1/(var_maxs-var_mins)
    C <- -var_mins/(var_maxs-var_mins)
    names(M) <- var_names
    names(C) <- var_names
    
    var_scales <- list(M = M, C = C)
    
    return(var_scales)
  }
  
  data2kdtree_WEP_2_parameter <- function(data_ddf, var_names, J_1, J_2, M, C){
    
    P = length(var_names)
    J = J_1*J_2
    K = (2*J+1)**P
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_names]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      scaled_data_subset <- t(M*t(data_subset)+C)
      T <- get_kdtree_summand_2_parameter(scaled_data_subset, J_1, J_2)
      for (i in 1:length(T)){collect(i, T[i])}
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- output + sum(unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_names = var_names, M = M, C = C, J_1 = J_1, J_2 = J_2)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys <- unlist(lapply(mr_out, function(x) x[[1]][1]))
    values <- unlist(lapply(mr_out, "[[", 2))
    
    kdtree_2D_power_SEP <- values[invPerm(as.integer(keys))]
    kdtree_2D_power_WEP <- kdtree_2D_power_SEP/kdtree_2D_power_SEP[1]
    
    kdtree_WEP <- array(kdtree_2D_power_WEP, dim=rep(2*J+1, P))
    
    f_tr <- function(WEP){
      res = Trig_2D_power_WEP2Trig_WEP(WEP, J_1, J_2)
      return(res)
    }
    
    for (p in 1:P){
      kdtree_WEP <- apply(kdtree_WEP, (1:(P-1)), f_tr)
    }
    
    return(kdtree_WEP)
  }
  
  data2kdtree_WEP_3_parameter <- function(data_ddf, var_names, J_1, J_2, J_3, M, C){
    
    P = length(var_names)
    J = J_1*J_2*J_3
    K = (2*J+1)**P
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i){return(data.frame(map.values[[i]][, var_names]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (length(data_subset) > 0){
      scaled_data_subset <- t(M*t(data_subset)+C)
      T <- get_kdtree_summand_3_parameter(scaled_data_subset, J_1, J_2, J_3)
      for (i in 1:length(T)){collect(i, T[i])}
    }
    })
    reduce <- expression(pre = {output <- 0}, reduce = {output <- output + sum(unlist(reduce.values))}, post = {collect(reduce.key, output)})
    parameter_list <- list(var_names = var_names, M = M, C = C, J_1 = J_1, J_2 = J_2, J_3=J_3)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    
    keys <- unlist(lapply(mr_out, function(x) x[[1]][1]))
    values <- unlist(lapply(mr_out, "[[", 2))
    
    kdtree_3D_power_SEP <- values[invPerm(as.integer(keys))]
    kdtree_3D_power_WEP <- kdtree_3D_power_SEP/kdtree_3D_power_SEP[1]
    
    kdtree_WEP <- array(kdtree_3D_power_WEP, dim=rep(2*J+1, P))
    
    f_tr <- function(WEP){
      res = Trig_3D_power_WEP2Trig_WEP(WEP, J_1, J_2, J_3)
      return(res)
    }
    
    for (p in 1:P){
      kdtree_WEP <- apply(kdtree_WEP, (1:(P-1)), f_tr)
    }
    
    return(kdtree_WEP)
  }
  
  Trig_2D_power_WEP2Trig_WEP <- function(Trig_2D_power_WEP, J_1, J_2) 
  {
    Trig_WEP_1 <- Trig_2D_power_WEP[1]
    Trig_2D_power_WEP <- matrix(Trig_2D_power_WEP[-1], nrow = 2)
    
    Cosine_2D_power_WEP <- array(Trig_2D_power_WEP[1, ], dim = c(J_1, J_2))
    Sine_2D_power_WEP <- array(Trig_2D_power_WEP[2, ], dim = c(J_1, J_2))
    
    Cosine_2D_product_WEP <- apply(apply(Cosine_2D_power_WEP, 1, ChevTrB), 1, ChevTrA)
    Sine_2D_product_WEP <- apply(apply(Sine_2D_power_WEP, 1, ChevTrB), 1, FUN = function(x) OddNegTr(ChevTrA(SineTr(x))))
    
    Cosine_WEP <- CosineSumTr(Cosine_2D_product_WEP)
    Sine_WEP <- SineSumTr(Sine_2D_product_WEP)
    
    Trig_WEP <- c(Trig_WEP_1, rbind(Cosine_WEP, Sine_WEP))
    
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
  
  WEP2multiplier_info <- function(WEP, lower, upper, d){
    P <- length(dim(WEP))
    d <- (d-1)%%P+1
    J <- as.integer((dim(WEP)[1]-1)/2)
    
    lower_sub <- lower[-d]
    upper_sub <- upper[-d]
    multiplier_array <- get_kdtree_multiplier(lower_sub[1], upper_sub[1], J)
    
    if (P > 2){
      for (p in 2:(P-1)){
        multiplier_array <- outer(multiplier_array, get_kdtree_multiplier(lower_sub[p], upper_sub[p], J))
      }
    }
    
    multiplier <- c()
    for (j in 1:(2*J+1)){
      multiplier <- c(multiplier, sum(asub(WEP, j, d)*multiplier_array))
    }
    multiplier_info <- list(d = d, multiplier = multiplier)
    
    return(multiplier_info)
  }
  
  WEP2median_info <- function(WEP, lower, upper, d){
    multiplier_info <- WEP2multiplier_info(WEP, lower, upper, d)
    d <- multiplier_info$d
    J <- as.integer((dim(WEP)[1]-1)/2)
    
    fn <- function(m){
      return(sum(multiplier_info$multiplier*(get_kdtree_multiplier(lower[multiplier_info$d], m, J)-get_kdtree_multiplier(m, upper[multiplier_info$d], J))))
    }
    
    median <- uniroot(f = fn, interval = c(lower[d], upper[d]))$root
    median_info <- list(d = d, median = median)
    return(median_info)
  }
  
  WEP2kdtree_info <- function(WEP, D, M, C){
    P <- length(dim(WEP))
    neighborhoods_unscaled <- list(list(lower = rep(0, P), upper = rep(1, P)))
    medians <- list()
    
    for(d in 1:D){
      neighborhoods_unscaled_info = WEP2sub_neighborhoods_unscaled_info(WEP = WEP, neighborhoods_unscaled = neighborhoods_unscaled, d = d, M = M, C = C)
      neighborhoods_unscaled = neighborhoods_unscaled_info$neighborhoods_unscaled
      medians[[d]] = neighborhoods_unscaled_info$medians
    }
    
    neighborhoods = scale_neighborhoods(neighborhoods_unscaled = neighborhoods_unscaled, M = M, C = C)
    kdtree <- list(neighborhoods = neighborhoods, medians = medians)
    
    return(kdtree)
  }
  
  WEP2sub_neighborhoods_unscaled_info <- function(WEP, neighborhoods_unscaled, d, M, C){
    new_neighborhoods_unscaled <- list()
    new_medians <- c()
    
    for(i in 1:length(neighborhoods_unscaled)){
      median_info <- WEP2median_info(WEP = WEP, lower = neighborhoods_unscaled[[i]]$lower, upper = neighborhoods_unscaled[[i]]$upper, d = d)
      d <- median_info$d
      new_medians <-c(new_medians, (median_info$median-C[d])/M[d])
      left_sub_neighborhood_lower <- neighborhoods_unscaled[[i]]$lower
      left_sub_neighborhood_upper <- neighborhoods_unscaled[[i]]$upper
      right_sub_neighborhood_lower <- neighborhoods_unscaled[[i]]$lower
      right_sub_neighborhood_upper <- neighborhoods_unscaled[[i]]$upper
      left_sub_neighborhood_upper[d] <- median_info$median
      right_sub_neighborhood_lower[d] <- median_info$median
      new_neighborhoods_unscaled[[2*i-1]] <- list(lower = left_sub_neighborhood_lower, upper = left_sub_neighborhood_upper)
      new_neighborhoods_unscaled[[2*i]] <- list(lower = right_sub_neighborhood_lower, upper = right_sub_neighborhood_upper)
    }
    
    neighborhoods_unscaled_info <- list(neighborhoods_unscaled = new_neighborhoods_unscaled, medians = new_medians)
    
    return(neighborhoods_unscaled_info)
  }
 
  scale_neighborhoods <- function(neighborhoods_unscaled, M, C){
    neighborhoods_scaled = neighborhoods_unscaled
    
    for(i in 1:length(neighborhoods_unscaled)){
      neighborhoods_scaled[[i]]$lower <- (neighborhoods_unscaled[[i]]$lower-C)/M
      neighborhoods_scaled[[i]]$upper <- (neighborhoods_unscaled[[i]]$upper-C)/M
    }
    
    return(neighborhoods_scaled)
  }
  
  kdtree_medians2cellcounts <- function (data_ddf, var_names, medians){
    D <- length(medians)
    medians_array <- matrix(0, D, 2^(D-1))
    
    for (d in 1:D){
      medians_array[d, 1:2^(d-1)] <- medians[[d]]
    }
    
    map <- expression({data_subset <- data.frame(data.table::rbindlist(lapply(seq_along(map.values), function(i) {return(data.frame(map.values[[i]][, var_names]))})))
    data_subset <- data_subset[complete.cases(data_subset), ]
    if (nrow(data_subset) > 0) {
      data_subset <- as.matrix(data_subset)
      C <- kdtree_medians2cellcounts(data_subset, medians_array)
      for (i in 1:length(C)){collect(i, C[i])}
    }
    })
    reduce <- expression(pre = {sum <- 0}, reduce = {sum <- sum + sum(unlist(reduce.values))}, post = {collect(reduce.key, sum)})
    parameter_list <- list(var_names = var_names, medians_array = medians_array)
    packages <- c("datadr", "data.table", "Rcpp", "drEPF")
    control <- rhipeControl(mapred = list(mapreduce.task.timeout = 0))
    
    mr_out <- getAttribute(mrExec(data_ddf, map = map, reduce = reduce, params = parameter_list, packages = packages, control = control), "conn")$data
    cell_counts <- unlist(lapply(mr_out, "[[", 2))
    return(cell_counts)
  }
  
  start_time <- proc.time()
  
  if(!inherits(data_ddf, "ddf")) stop("data_ddf should be a distributed data frame.")
  
  for (var_name in var_names){
    if(!is.character(var_name) | var_name %in% names(attributes(data_ddf)$ddf$vars) == FALSE | attributes(data_ddf)$ddf$vars[var_name] != "numeric") stop("var_names should contain characters representing the numeric-variables present in the data to construct the KD-tree.")
  }
  
  if(!is.double(D)) stop("D should be an integer representing the depth of the KD-tree.")
  
  if(!D == floor(D)) stop("D should be an integer representing the depth of the KD-tree.")
  
  if(is.null(var_ranges) == TRUE) {
    var_scales = data2var_scales(data_ddf, var_names)
    M = var_scales$M
    C = var_scales$C
  } else if(!is.vector(var_ranges[[1]]) | !is.double(var_ranges[[1]]) | !length(var_ranges[[1]]) == length(var_names) | !is.vector(var_ranges[[2]]) | !is.double(var_ranges[[2]]) | !length(var_ranges[[2]]) == length(var_names)) {
    stop("var_ranges should contain the minimums and maximuams of the variable in var_names")
  } else if(!all(var_ranges[[1]] < var_ranges[[2]])) {
    stop("var_ranges should contain the minimums and maximuams of the variable in var_names")
  } else {
    M = 1/(var_ranges[[2]]-var_ranges[[1]])
    C = -var_ranges[[1]]/(var_ranges[[2]]-var_ranges[[1]])
  }
  
  if(!is.vector(J) | !all(J==floor(J)) | length(J) %in% c(2, 3) == FALSE) stop("J should be and integer vector containing 3 or 4 parameters.")
  
  if(length(J) == 2) kdtree_WEP = data2kdtree_WEP_2_parameter(data_ddf, var_names, J_1 = J[1], J_2 = J[2], M, C)
  
  if(length(J) == 3) kdtree_WEP = data2kdtree_WEP_3_parameter(data_ddf, var_names, J_1 = J[1], J_2 = J[2], J_3 = J[3], M, C)
  
  kdtree_info <- WEP2kdtree_info(kdtree_WEP, D, M, C)
  
  tree <- list()
  for (d in 1:D){
    tree[[d]] = list(Level = d, Splitting_axis = names(kdtree_info$medians[[d]])[1], Medians = unname(kdtree_info$medians[[d]]))
  }
  
  if(include_cell_counts == TRUE){
    counts = kdtree_medians2cellcounts(data_ddf, var_names, kdtree_info$medians)
  }
  
  stop_time = proc.time()
  
  if(include_cell_counts == TRUE){
    kdtree = list(Tree = tree, Neighborhoods = kdtree_info$neighborhoods, Counts = counts, Computation_time_in_seconds = unname(stop_time[3] - start_time[3]))
  } else {
    kdtree = list(Tree = tree, Neighborhoods = kdtree_info$neighborhoods, Computation_time_in_seconds = unname(stop_time[3] - start_time[3]))
  }
  
  return(kdtree)
}