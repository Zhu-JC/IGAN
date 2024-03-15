#' Grouping spots pairs into groups
#'
#' @param spatial spatial coordinate information, a dataframe which rows represent spots and cols represent coordinate
#' @param ident celltype information, a dataframe which rows represent spots and cols represent celltype
#' @param group_size The group size parameters determined by humans. For single-cell spatial RNA-seq data, we recommend grouping them into 5000 spots, while for spatial bulk RNA-seq data, grouping into 2000 spots is recommended.
#'
#' @return grouping results in spot pairs, each group of 5 rows represents:
#'         spatial distance between spot pairs, index of the spot sending the signal, index of the spot receiving the signal, cell type of the spot sending the signal, cell type of the spot receiving the signal.
#'
#' @export
#'
#' @examples cells_group = cells_to_group(spatial, ident, group_size)
cells_to_group = function(spatial, ident, group_size, celltype = NULL){
  ncells = nrow(spatial)
  near_dis = c()
  for (i in 1:ncells) {
    dis = sqrt((spatial[,1] - spatial[i,1])^2 + (spatial[,2] - spatial[i,2])^2)
    dis = dis[which(dis>0)]
    near_dis[i] = min(dis)
  }

  rank_dis = sort(near_dis)
  threshold_dis = rank_dis[round(0.99*ncells)]
  if(is.null(celltype)){
    celltype = unique(ident)
  }

  all = apply(celltype, 1, function(x, spatial, ident, threshold_dis){
    cell1 = which(ident[,1] == x[1])
    all1 = apply(celltype, 1, function(y, cell1 , spatial, ident, threshold_dis){
      cell2 = which(ident[,1] == y[1])
      spatial1 = as.matrix(spatial[cell1,])
      spatial2 = as.matrix(spatial[cell2,])
      all2 = apply(spatial1, 1, function(z, spatial2, threshold_dis, cell1, cell2, x, y, spatial1){
        dis = sqrt((spatial2[,1] - z[1])^2 + (spatial2[,2] - z[2])^2)
        use_dis = which(dis>0 & dis<1.25*threshold_dis)
        all3 = as.data.frame(matrix(nrow = 5, ncol = length(use_dis)))
        all3[1,] = dis[use_dis]
        all3[2,] = cell1[which(spatial1[,1] == z[1] & spatial1[,2] == z[2])]
        all3[3,] = cell2[use_dis]
        all3[4,] = x[1]
        all3[5,] = y[1]
        return(all3)
      }, spatial2, threshold_dis, cell1, cell2, x, y, spatial1)
      return(as.data.frame(all2))

    }, cell1, spatial, ident, threshold_dis)
    return(all1)

  }, spatial, ident, threshold_dis)

  all2 = list()
  for (i in 1:length(all)) {
    all2 = c(all2, all[[i]])

  }
  names(all2) = 1:length(all2)
  sz = as.data.frame(sapply(all2, function(x) {
    return(ncol(x))
  }))
  all2 = all2[-which(sz == 0)]

  all3 = list()
  s = 1
  u = 0
  l_all2 = length(all2)
  for (i in 1:l_all2) {
    if(length(all2) > 0){
      while (ncol(all2[[i-u]]) > 1.1*group_size) {
        all3[[s]] = all2[[i-u]][,1:group_size]
        all2[[i-u]] = all2[[i-u]][,-(1:group_size)]
        s = s+1
      }

      if(ncol(all2[[i-u]]) <= 1.1*group_size & ncol(all2[[i-u]]) >= 0.9*group_size){
        all3[[s]] = all2[[i-u]]
        all2 = all2[-(i-u)]
        s = s+1
        u = u+1
      }else{
        sz = as.data.frame(sapply(all2, function(x) {
          return(ncol(x))
        }))
        if(nrow(sz)>1){
          min_sz = which(sz == min(sz[-1,1]))
          all2[[i-u]] = cbind(all2[[i-u]], all2[i-u-1+min_sz])
          all2 = all2[-(i-u-1+min_sz)]
        }
        u = u+1
      }
    }
  }
  if(length(all2) > 0){
    if(ncol(all2[[1]])<0.6*group_size){
      sz = as.data.frame(sapply(all3, function(x) {
        return(ncol(x))
      }))
      min_sz = which(sz == min(sz))
      all3[[min_sz]] = cbind(all3[[min_sz]],all2[[1]])
    }else{
      all3[[s]] = all2[[1]]
    }
  }
  return(all3)
}















