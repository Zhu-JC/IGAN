#' Get the cellchat's network circle graph
#'
#' @param useGene the gene's name used to plot
#' @param result result get from function IGAN
#' @param cells_group grouping result from function cells_to_group or cells_to_group10X
#' @param gene gene name of the dataset
#' @param match_list1 gene list selected in sending spots
#' @param ident celltype information, a dataframe which rows represent spots and cols represent celltype
#'
#' @return a cellchat's network circle graph
#' @export
#'
#' @examples plot_cellchat('SPP1',result,cells_group,gene,match_list1,ident)
plot_cellchat = function(useGene,result,cells_group,gene,match_list1,ident){
ncelltype = length(unique(ident[,1]))
celltypes = unique(ident[,1])
cpsize = matrix(0,nrow = ncelltype,ncol = ncelltype)
ngroup = length(cells_group)
for (i in 1:ngroup){
  sub_group = cells_group[[i]]
  for (j in 1:ncol(sub_group)) {
    cp1 = which(celltypes == sub_group[4,j])
    cp2 = which(celltypes == sub_group[5,j])
    cpsize[cp1,cp2] = cpsize[cp1,cp2] + 1
  }
}
cpsize = cpsize+1

gene = toupper(gene)
useGene = toupper(useGene)
use_rank = which(match_list1 == match(useGene, gene))

celltypes_corM = matrix(0,nrow = ncelltype,ncol = ncelltype)
row.names(celltypes_corM) = celltypes
colnames(celltypes_corM) = celltypes
for (i in 1:ngroup) {
  use_celltype = as.matrix(cells_group[[i]])[4:5,]
  sub_result = result[[i]][use_rank]
  k = which(colSums(sub_result[[1]])>0)
  cell = use_celltype[,as.numeric(colnames(sub_result[[1]])[k])]
  for (j in 1:ncol(cell)) {
    c1 = which(celltypes == cell[1,j])
    c2 = which(celltypes == cell[2,j])
    celltypes_corM[c1,c2] = celltypes_corM[c1,c2] + 1
  }
}
celltypes_corM = celltypes_corM/cpsize
celltype_size = table(ident)[match(celltypes,names(table(ident)))]
cellchat_graph = netVisual_circle(celltypes_corM, vertex.weight = celltype_size, weight.scale = T, label.edge= F, title.name = useGene)
return(cellchat_graph)
}


