
#' Get the graph of CCC pattern
#'
#' @param result result get from function IGAN
#' @param cells_group grouping result from function cells_to_group or cells_to_group10X
#' @param gene gene name of the dataset
#' @param match_list1 gene list selected in sending spots
#' @param match_list2 gene list selected in receiving spots
#' @param ident celltype information, a dataframe which rows represent spots and cols represent celltype
#' @param pattern set it's sending pattern or receiving pattern
#' @param n the number of pattern
#'
#' @return graph of CCC pattern
#' @export
#'
#' @examples pattern_graph = plot_pattern(result,cells_group,gene,match_list1,match_list2,ident,'outgoing',4)
plot_pattern = function(result,cells_group,gene,match_list1,match_list2,ident,pattern,n){
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

  #get_cpsize
  ncelltype = length(unique(ident[,1]))
  celltypes = unique(ident[,1])
  cpsize = matrix(0,nrow = ncelltype,ncol = ncelltype)
  for (i in 1:ngroup){
    sub_group = cells_group[[i]]
    for (j in 1:ncol(sub_group)) {
      cp1 = which(celltypes == sub_group[4,j])
      cp2 = which(celltypes == sub_group[5,j])
      cpsize[cp1,cp2] = cpsize[cp1,cp2] + 1
    }
  }
  cpsize = cpsize+1

  celltype_gene = matrix(0,nrow = ncelltype,ncol = length(use_gene))
  for (i in 1:ngroup) {
    #print(i)
    sub_result = result[[i]][use_gene]
    if(pattern == 'outgoing'){
      use_celltype = as.matrix(cells_group[[i]])[4:5,]
    }else{
      use_celltype = as.matrix(cells_group[[i]])[5:4,]
    }
    group_cor = sapply(sub_result[lengths(sub_result) > 0], function(x, use_celltype, sz, celltypes){
      k = which(colSums(x)>0)
      cell = use_celltype[,as.numeric(colnames(x)[k])]
      celltypes_corM = matrix(0,nrow = sz,ncol = sz)
      for (j in 1:ncol(cell)) {
        c1 = which(celltypes == cell[1,j])
        c2 = which(celltypes == cell[2,j])
        celltypes_corM[c1,c2] = celltypes_corM[c1,c2] + 1
      }
      celltypes_corM = celltypes_corM/cpsize
      return(rowSums(celltypes_corM))
    }, use_celltype, sz = ncelltype, celltypes = unique(ident[,1]))
    celltype_gene[,lengths(sub_result) > 0] = celltype_gene[,lengths(sub_result) > 0] + group_cor
  }

  celltype_gene = as.data.frame(celltype_gene)
  colnames(celltype_gene) = names(use_gene)
  row.names(celltype_gene) = celltypes
  celltype_gene = celltype_gene[,colSums(celltype_gene)!=0]
  celltype_gene = sweep(celltype_gene, 2L, apply(celltype_gene, 2, function(x) max(x, na.rm = TRUE)), '/', check.margin = FALSE)
  data = as.matrix(celltype_gene)
  options(warn = -1)
  data = data[rowSums(data)!=0,]
  title.name = paste0(pattern, " signaling \n")

  do.facet = TRUE
  k.range = seq(2,10)
  nrun = 30
  seed.use = 10

  res <- NMF::nmfEstimateRank(data, range = k.range, method = 'lee', nrun=nrun, seed = seed.use)
  df1 <- data.frame(k = res$measures$rank, score = res$measures$cophenetic, Measure = "Cophenetic")
  df2 <- data.frame(k = res$measures$rank, score = res$measures$silhouette.consensus, Measure = "Silhouette")
  df <- rbind(df1, df2)

  k = n
  color.use = NULL
  color.heatmap = "Spectral"
  title.legend = "Contributions"
  width = 6
  height = 15
  font.size = 8

  outs_NMF <- NMF::nmf(data, rank = k, method = 'lee', seed = 'nndsvd')
  W <- scaleMat(outs_NMF@fit@W, 'r1')
  H <- scaleMat(outs_NMF@fit@H, 'c1')
  colnames(W) <- paste0("Pattern ", seq(1,ncol(W)));
  rownames(H) <- paste0("Pattern ", seq(1,nrow(H)));

  net <- W
  if (is.null(color.use)) {
    color.use <- scPalette(length(rownames(net)))
  }
  color.heatmap = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(255)

  df<- data.frame(group = rownames(net)); rownames(df) <- rownames(net)
  cell.cols.assigned <- setNames(color.use, unique(as.character(df$group)))
  row_annotation <- HeatmapAnnotation(df = df, col = list(group = cell.cols.assigned),which = "row",
                                      show_legend = FALSE, show_annotation_name = FALSE,
                                      simple_anno_size = grid::unit(0.2, "cm"))

  ht1 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                left_annotation = row_annotation,
                cluster_rows = F,cluster_columns = F,clustering_method_rows = "average",
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                show_heatmap_legend = F,
                column_title = "Cell patterns",column_title_gp = gpar(fontsize = 15)
  )


  net <- t(H)

  ht2 = Heatmap(net, col = color.heatmap, na_col = "white", name = "Contribution",
                cluster_rows = T,cluster_columns = F,clustering_method_rows = "average",
                row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                width = unit(width, "cm"), height = unit(height, "cm"),
                column_title = "Communication patterns",column_title_gp = gpar(fontsize = 15),
                heatmap_legend_param = list(title = title.legend, title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                            border = NA, at = c(round(min(net, na.rm = T), digits = 1), round(max(net, na.rm = T), digits = 1)),
                                            legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 6),grid_width = unit(2, "mm"))
  )

  gb_ht1 = grid.grabExpr(draw(ht1))
  gb_ht2 = grid.grabExpr(draw(ht2))

  pattern_graph = list(gb_ht1,gb_ht2)
  return(pattern_graph)

}


