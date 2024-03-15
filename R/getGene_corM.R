#' Title
#'
#' @param result
#' @param cells_group
#' @param gene
#' @param match_list1
#' @param match_list2
#'
#' @return
#' @export
#'
#' @examples
getGene_corM = function(result, cells_group, gene, match_list1, match_list2){
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
}
