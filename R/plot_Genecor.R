
#' Get the graph of useGene's CCC activity of every spot
#'
#' @param useGene the gene's name used to plot
#' @param result result get from function IGAN
#' @param cells_group grouping result from function cells_to_group or cells_to_group10X
#' @param spatial spatial coordinate information, a dataframe which rows represent spots and cols represent coordinate
#' @param match_list1 gene list selected in sending spots
#' @param gene gene name of the dataset
#'
#' @return graph of useGene's CCC activity of every spot
#' @export
#'
#' @examples a = plot_Genecor('SPP1', result, cells_group, spatial, match_list1, gene)
#'           a + scale_color_distiller(palette = "Set1")

plot_Genecor = function(useGene, result, cells_group, spatial, match_list1, gene){
  deep = rowSums(sapply(cells_group, function(x,sz){
    group_member = rep(0, sz)
    a = as.data.frame(table(as.matrix(x[2,])))
    group_member[as.numeric(levels(a[,1]))] = group_member[as.numeric(levels(a[,1]))] + a[,2]
    return(group_member)
  }, sz = nrow(spatial)))
  deep = deep+1

  gene = toupper(gene)
  useGene = toupper(useGene)
  use_rank = which(match_list1 == match(useGene, gene))

  ngroup = length(cells_group)
  all_cor = rep(0,nrow(spatial))
  for (i in 1:ngroup) {
    #print(i)
    sub_result = result[[i]][use_rank]
    use_cell = as.matrix(cells_group[[i]])[2:3,]
    group_cor = rowSums(sapply(sub_result, function(x, use_cell, sz){
      k = which(colSums(x)>0)
      cells_sent = as.numeric(use_cell[1, as.numeric(colnames(x)[k])])
      cells_sent = as.matrix(as.data.frame(table(cells_sent)))
      cells_rec = as.numeric(use_cell[2,as.numeric(colnames(x)[k])])
      cells_rec = as.matrix(as.data.frame(table(cells_rec)))
      sub_cor = rep(0,sz)
      sub_cor[as.numeric(cells_sent[,1])] = sub_cor[as.numeric(cells_sent[,1])] + as.numeric(cells_sent[,2])
      sub_cor[as.numeric(cells_rec[,1])] = sub_cor[as.numeric(cells_rec[,1])] + as.numeric(cells_rec[,2])
      return(sub_cor)
    }, use_cell, sz = nrow(spatial)))
    all_cor = all_cor+group_cor
  }
  all_cor = all_cor/deep
  all_cor = log(all_cor+1)
  data <- data.frame(x = spatial[,1],
                     y = spatial[,2],
                     color = all_cor)
  G = ggplot(data, aes(x, y, color = color)) +
    geom_point() +
    labs(title = paste0(useGene,' CCC'), x = "X", y = "Y")
  return(G)
}

