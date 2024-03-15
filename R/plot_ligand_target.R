#' Title
#'
#' @param result
#' @param spatial
#' @param cells_group
#' @param gene
#' @param match_list1
#' @param match_list2
#' @param genename
#'
#' @return
#' @export
#'
#' @examples
plot_ligand_target = function(result, spatial, cells_group, gene, match_list1, match_list2, genename){
  useGene = genename
  ngroup = length(cells_group)

  data1 <- data.frame(x = spatial[,1],
                      y = spatial[,2])
  data1[,3:5] = 0
  data2 <- data.frame(x = spatial[,1],
                      y = spatial[,2])
  data2[,3:4] = 0
  use_gene = match(useGene, gene[match_list1])
  for (j in 1:ngroup) {
    #print(i)
    sub_result = result[[j]][[use_gene]]
    if(is.null(sub_result[[1]])){
      next
    }
    target = as.integer(row.names(sub_result))
    for (v in 1:nrow(sub_result)) {
      k = which(sub_result[v,]>0)
      k = as.integer(colnames(sub_result)[k])
      sender = as.integer(cells_group[[j]][2,k])
      receiver = as.integer(cells_group[[j]][3,k])

      data1[sender,3] = data1[sender,3] + t(M[match(useGene,gene),sender])
      data1[sender,4] = data1[sender,4] + t(M[match_list2[target[v]],receiver])
      data1[sender,5] = data1[sender,5] + 1
      data2[receiver,3] = data2[receiver,3] + t(M[match_list2[target[v]],receiver])
      data2[receiver,4] = data2[receiver,4] + 1
    }


  }
  use = which(data1[,5]>0)
  data1[use,3] = data1[use,3]/data1[use,5]
  data1[use,4] = data1[use,4]/data1[use,5]

  r = cor(data1[use,3],data1[use,4],method = 'spearman')


  colnames(data1) = c('x', 'y', 'color', 'neighbor')
  grey = which(data1[,5] == 0)

  data_gray <- data1[grey,]
  data_gradient <- data1[use,]

  G1 = ggplot() +
    geom_point(data = data_gray, aes(x = x, y = y), color = "gray") +
    geom_point(data = data_gradient, aes(x = x, y = y, color = color)) +
    scale_color_gradient(low = "lightpink", high = "darkred") +
    theme_minimal() +
    labs(title = paste0(useGene,' expression'), x = "X", y = "Y") +
    theme(text = element_text(size = 20),
          plot.title = element_text(size = 22),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))
  #G1 = ggplot(data1) +
    #geom_point(aes(x = x, y = y, color = color), data = data1[grey,]) +
    #geom_point(aes(x = x, y = y, color = color), data = data1[use,]) +
    #scale_color_gradient(low = "lightpink", high = "darkred", guide = 'none') +
    #scale_color_manual(values = 'gray', guide = 'none', aesthetics = 'color') +
    #theme_minimal() +
    #labs(title = paste0(useGene,' expression'), x = "X", y = "Y")+
    #theme(text = element_text(size = 20),  # 调整全局文本大小
          #plot.title = element_text(size = 22), # 调整标题大小
          #axis.title.x = element_text(size = 18), # 调整X轴标题大小
          #axis.title.y = element_text(size = 18)) # 调整Y轴标题大小
    #guides(color = guide_legend(title = "Intensity"))


  colnames(data2) = c('x', 'y', 'color')
  use2 = which(data2[,4]>0)
  data2[use2,3] = data2[use2,3]/data2[use2,4]
  grey = which(data2[,4] == 0)

  data_gray <- data2[grey,]
  data_gradient <- data2[use2,]

  G2 = ggplot() +
    geom_point(data = data_gray, aes(x = x, y = y), color = "gray") +
    geom_point(data = data_gradient, aes(x = x, y = y, color = color)) +
    scale_color_gradient(low = "lightpink", high = "darkred") +
    theme_minimal() +
    labs(title = paste0(useGene,' target genes expression'), x = "X", y = "Y") +
    theme(text = element_text(size = 20),
          plot.title = element_text(size = 22),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18))

  #G2 = ggplot(data2, aes(x, y, color = color)) +
    #geom_point() +
    #labs(title = paste0(useGene,' target genes expression'), x = "X", y = "Y")+
    #theme(text = element_text(size = 20),  # 调整全局文本大小
          #plot.title = element_text(size = 22), # 调整标题大小
          #axis.title.x = element_text(size = 18), # 调整X轴标题大小
          #axis.title.y = element_text(size = 18)) # 调整Y轴标题大小
  r = round(r,2)

  graph = grid.arrange(G1, G2, ncol = 2, top = textGrob(paste0('correlation = ',r), gp = gpar(fontsize = 20)))
  return(graph)
}
