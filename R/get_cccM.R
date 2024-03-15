
#' Get the matrix contain every spots' CCC feature
#'
#' @param result result get from function IGAN
#' @param cells_group grouping result from function cells_to_group or cells_to_group10X
#' @param spatial spatial coordinate information, a dataframe which rows represent spots and cols represent coordinate
#' @param gene gene name of the dataset
#' @param match_list1 gene list selected in sending spots
#' @param match_list2 gene list selected in receiving spots
#'
#' @return a matrix contain every spots' CCC feature
#' @export
#'
#' @examples cccM = get_cccM(result, cells_group, spatial, gene, match_list1, match_list2)
get_cccM = function(result, cells_group, spatial, gene, match_list1, match_list2){
  deep = rowSums(sapply(cells_group, function(x,sz){
    group_member = rep(0, sz)
    a = as.data.frame(table(as.matrix(x[2,])))
    group_member[as.numeric(levels(a[,1]))] = group_member[as.numeric(levels(a[,1]))] + a[,2]
    return(group_member)
  }, sz = nrow(spatial)))
  deep = deep+1

  ngroup = length(cells_group)
  Gene_corM = matrix(0,nrow = length(match_list1), ncol = length(match_list2))
  rownames(Gene_corM) = gene[match_list1]
  colnames(Gene_corM) = gene[match_list2]
  for (i in 1:ngroup) {
    #print(i)
    sub_result = result[[i]][lengths(result[[i]]) > 0]
    group_corM = sapply(sub_result, function(x, ncols){
      k = as.numeric(row.names(x))
      sub_corM = matrix(0,1,ncols)
      sub_corM[1,k] = 1
      return(sub_corM)
    }, ncols = length(match_list2))
    group_corM = t(group_corM)
    Gene_corM[lengths(result[[i]]) > 0,] = Gene_corM[lengths(result[[i]]) > 0,] + group_corM
  }
  use_gene = which(rowSums(Gene_corM>0) > 0)

  CCCM = matrix(0,nrow = nrow(spatial),ncol = 2*length(use_gene))
  for (i in 1:ngroup) {
    #print(i)
    sub_result = result[[i]][use_gene]
    use_cell = as.matrix(cells_group[[i]])[2:3,]
    group_M = sapply(sub_result[lengths(sub_result) > 0], function(x, use_cell, sz){
      k = which(colSums(x)>0)
      sumX = colSums(x)[k]
      cells_sent = as.numeric(use_cell[1, as.numeric(colnames(x)[k])])
      sent_cor = rep(0,sz)
      cells_rec = as.numeric(use_cell[2,as.numeric(colnames(x)[k])])
      rec_cor = rep(0,sz)
      for (j in 1:length(k)) {
        sent_cor[cells_sent[j]] = sent_cor[cells_sent[j]] + sumX[j]
        rec_cor[cells_rec[j]] = rec_cor[cells_rec[j]] + sumX[j]
      }

      return(cbind(sent_cor,rec_cor))
    }, use_cell, sz = nrow(spatial))
    CCCM[,lengths(sub_result) > 0] = group_M[1:nrow(spatial),]
    CCCM[,which((lengths(sub_result)>0) == 'TRUE')+length(use_gene)] = group_M[(nrow(spatial)+1):(2*nrow(spatial)),]
  }
  CCCM = CCCM/deep
  CCCM = log(CCCM+1)
  return(CCCM)
}
