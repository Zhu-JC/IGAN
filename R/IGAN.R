
#' Compute gene-gene associations in every single spots pair
#'
#' @param match_list1 gene list selected in sending spots
#' @param match_list2 gene list selected in receiving spots
#' @param gene gene name of the dataset
#' @param cells_group grouping result from function cells_to_group or cells_to_group10X
#' @param M gene expression matrix, a dataframe which cols represent genes and rows represent spots
#' @param p_value p_value used in Statistical Testing
#'
#' @return
#' @export
#'
#' @examples result = IGAN(match_list1, match_list2, gene, cells_group, M, p_value)
IGAN = function(match_list1, match_list2, gene, cells_group, M, p_value){

  CCC = sapply(cells_group, function(x, match_list1, match_list2, M, p_value){

    n_cp = ncol(x)
    M1 = M[match_list1,as.numeric(x[2,])]
    use_row = which(rowSums(M1>0) > 20)
    M1 = M1[use_row, ]
    M2 = M[match_list2,as.numeric(x[3,])]
    size_M1 = dim(M1)
    size_M2 = dim(M2)

    start_time <- Sys.time()
    neighbor1 = apply(M1, 1, function(x, size_M1){
      s1 = sort(x)
      s2 = match(s1,x)
      n0 = size_M1[2] - sum(sign(s1))
      h = round(0.05*sum(sign(s1)))
      k = n0+1
      upper = matrix(0, nrow = 1, ncol = size_M1[2])
      lower = matrix(0, nrow = 1, ncol = size_M1[2])
      while (k <= size_M1[2]) {
        #s = 0
        #while (k+s+1 <= size_M1[2] & s1[k+s+1] == s1[k]) {
        #  s = s+1
        #}
        neighbors = which(x == s1[k])
        s = length(which(s1 == s1[k]))-1
        if(s >= h){
          upper[1,neighbors] = x[s2[k]]
          lower[1,neighbors] = x[s2[k]]
        }else{
          upper[1,neighbors] = x[s2[min(size_M1[2],k+s+h)]]
          lower[1,neighbors] = x[s2[max(n0*(n0>h)+1, k-h)]]
        }
        k = k+s+1
      }
      neighbor1 = list(upper, lower)
      return(neighbor1)
    }, size_M1)
    end_time1 <- Sys.time()
    print(end_time1 - start_time)

    start_time <- Sys.time()
    neighbor2 = apply(M2, 1, function(x, size_M2){
      s1 = sort(x)
      s2 = match(s1,x)
      n0 = size_M2[2] - sum(sign(s1))
      h = round(0.05*sum(sign(s1)))
      k = n0+1
      upper = matrix(0, nrow = 1, ncol = size_M2[2])
      lower = matrix(0, nrow = 1, ncol = size_M2[2])
      while (k <= size_M2[2]) {
        #s = 0
        #while (k+s+1 <= size_M1[2] & s1[k+s+1] == s1[k]) {
        #  s = s+1
        #}
        neighbors = which(x == s1[k])
        s = length(which(s1 == s1[k]))-1
        if(s >= h){
          upper[1,neighbors] = x[s2[k]]
          lower[1,neighbors] = x[s2[k]]
        }else{
          upper[1,neighbors] = x[s2[min(size_M1[2],k+s+h)]]
          lower[1,neighbors] = x[s2[max(n0*(n0>h)+1, k-h)]]
        }
        k = k+s+1
      }
      neighbor2 = list(upper, lower)
      return(neighbor2)
    }, size_M2)
    end_time1 <- Sys.time()
    print(end_time1 - start_time)

    upper1 = do.call(rbind, lapply(neighbor1, function(x) x[[1]]))
    lower1 = do.call(rbind, lapply(neighbor1, function(x) x[[2]]))
    rm(neighbor1)
    upper2 = do.call(rbind, lapply(neighbor2, function(x) x[[1]]))
    lower2 = do.call(rbind, lapply(neighbor2, function(x) x[[2]]))
    rm(neighbor2)

    p = qnorm(1-p_value, 0, 1)
    n_list1 = matrix(1, nrow = 1, ncol = length(match_list1))
    DeGe = lapply(1:length(match_list1), function(i) data.frame())
    #DeCP = lapply(1:length(match_list1), function(i) data.frame())

    start_time <- Sys.time()
    for (i in 1:size_M2[2]) {
      g1 = which(M1[,i] >0)
      g2 = which(M2[,i] >0)

      lower1_g1i = lower1[g1, i]
      upper1_g1i = upper1[g1, i]
      lower2_g2i = lower2[g2, i]
      upper2_g2i = upper2[g2, i]
      neighbor_x = +(M1[g1, ] <= matrix(upper1_g1i, nrow = length(g1), ncol = size_M2[2]) & M1[g1, ] >= matrix(lower1_g1i, nrow = length(g1), ncol = size_M2[2]))
      neighbor_y = +(M2[g2, ] <= matrix(upper2_g2i, nrow = length(g2), ncol = size_M2[2]) & M2[g2, ] >= matrix(lower2_g2i, nrow = length(g2), ncol = size_M2[2]))

      n_x = rowSums(neighbor_x)
      n_y = rowSums(neighbor_y)
      cp_network = (neighbor_x %*% t(neighbor_y) * n_cp - n_x %*% t(n_y))/sqrt((n_x %*% t(n_y))*((n_x - n_cp) %*% t(n_y - n_cp))/(n_cp-1))

      cp_network = ifelse((cp_network>p),1,0)
      corr_g = which(rowSums(cp_network) > 0)
      for (j in corr_g) {
        DG = which(cp_network[j,] > 0)
        DeGe[[use_row[g1[j]]]][1:length(DG), n_list1[use_row[g1[j]]]] = g2[DG]
        colnames(DeGe[[use_row[g1[j]]]])[n_list1[use_row[g1[j]]]] = i
        #DeCP[[use_row[g1[j]]]][1, n_list1[use_row[g1[j]]]] = i
        n_list1[use_row[g1[j]]] = n_list1[use_row[g1[j]]] + 1
      }
    }
    DeGe2 = sapply(DeGe, function(x, size_M2){
      if(nrow(x)>0){
        target_M = apply(x, 2, function(x,size_M2){
          gg = matrix(0,1,size_M2[1])
          gg[x[which(x>0)]] = 1
          return(t(gg))
        }, size_M2)
        target_M = as.data.frame(target_M)
        target = which(rowSums(target_M) > 0.01*size_M2[2])
        if(length(target)>0){
          target_M2 = target_M[target,]
          return(target_M2)
        }
      }
    }, size_M2)

    result = list(DeGe2)
    end_time <- Sys.time()
    print(end_time - start_time)
    return(result)
  }, match_list1, match_list2, M, p_value)

}
