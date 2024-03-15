#' Get the  sankey graph of ligand's GO-ligand-recepter-downstream pathway
#'
#' @param result result get from function IGAN
#' @param cells_group grouping result from function cells_to_group or cells_to_group10X
#' @param gene gene name of the dataset
#' @param match_list1 gene list selected in sending spots
#' @param match_list2 gene list selected in receiving spots
#' @param p p_value used in KEGG enriched analysis
#' @param minGenes the min number of associated genes in ligands' downstream pathways
#'
#' @return a sankey graph of ligand's GO-ligand-recepter-downstream pathway
#' @export
#'
#' @examples sankey_graph = plot_sankey(result,cells_group,gene,match_list1,match_list2,0.001,5)
plot_sankey = function(result,cells_group,gene,match_list1,match_list2,p,minGenes, OrgDb, send_cell = NULL, rec_cell = NULL, threshold){
  ngroup = length(result)
  result2 = vector("list", length = ngroup)
  if(!is.null(send_cell) & !is.null(rec_cell)){
    for (i in 1:ngroup) {
      usel = which(lengths(result[[i]]) > 0)
      usep = which(cells_group[[i]][4,] == send_cell & cells_group[[i]][5,] == rec_cell)
      sub_result = vector("list", length = length(result[[i]]))
      for (j in 1:length(usel)) {
        subp = na.omit(match(usep,colnames(result[[i]][[usel[j]]])))
        sub_result[[usel[j]]] = result[[i]][[usel[j]]][,subp]
      }
      result2[[i]] = sub_result
    }
  }else{
    result2 = result
  }

  Gene_corM = matrix(0,nrow = length(match_list1), ncol = length(match_list2))
  rownames(Gene_corM) = gene[match_list1]
  colnames(Gene_corM) = gene[match_list2]
  for (i in 1:ngroup) {
    #print(i)
    sub_result = result2[[i]][lengths(result2[[i]]) > 0]
    if(length(sub_result)>0){
      group_corM = sapply(sub_result, function(x, ncols){
        k = as.numeric(row.names(x))
        sub_corM = matrix(0,1,ncols)
        sub_corM[1,k] = 1
        return(sub_corM)
      }, ncols = length(match_list2))
      group_corM = t(group_corM)
      Gene_corM[lengths(result2[[i]]) > 0,] = Gene_corM[lengths(result2[[i]]) > 0,] + group_corM
    }
  }
  use_gene = which(rowSums(Gene_corM>0) > threshold)

  ligand = as.data.frame(gene[match_list1[use_gene]])
  ligand_id <- bitr(ligand[,1], fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = OrgDb)

  ego <- enrichGO(gene          = ligand_id$ENTREZID,
                  OrgDb         = OrgDb,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  use_go = which(ego@result[["p.adjust"]] < p)
  go_result = ego@result[["geneID"]][use_go]
  go_result = strsplit(go_result,'/')
  go_name = ego@result[["Description"]][use_go]

  #开始准备桑基图
  if(metadata(OrgDb)[6,2] == "Human"){
    CellChatDB <- CellChatDB.human
    org = 'hsa'
  }else{
    CellChatDB <- CellChatDB.mouse
    org = 'mmu'
  }
  useligand = ligand[,1]
  shangyou = data.frame(matrix(ncol = 1))
  xiayou = data.frame(matrix(ncol = 1))
  R = c()
  colnames(shangyou) = 'name'
  h = 1
  p = 1
  for (v in 1:length(useligand)) {
    data = colnames(Gene_corM)[which(Gene_corM[use_gene[v],]>0)]
    nodes = c()
    s = 1
    for (i in 1:length(use_go)) {
      if(useligand[v] %in% go_result[[i]]){
        nodes[s] = go_name[i]
        s = s+1
      }

    }
    if(s == 1){
      next
    }
    k = s-1
    nodes = as.data.frame(nodes)
    nodes['ligand',1] = useligand[v]
    colnames(nodes) = 'name'

    links = data.frame('source' = c(0:(k-1)), 'target' = rep(k,k), 'value' = rep(0,k))
    color = letters
    # color = c("red", "green", "blue")
    data = as.character(as.matrix(data))
    test = bitr(data, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
    kk <- enrichKEGG(gene = test$ENTREZID,
                     organism = org,
                     pvalueCutoff = 1,
                     qvalueCutoff = 1)
    if(is.null(kk)){
      next
    }
    #找受体
    rank = which(CellChatDB[["interaction"]][["ligand"]] == useligand[v])
    recepter_name = CellChatDB[["interaction"]][["receptor.symbol"]][rank]
    recepter_name = recepter_name[!grepl(':',recepter_name)]
    recepter = strsplit(recepter_name,', ')
    allSymbols <- keys(OrgDb, keytype="SYMBOL")
    recepter = sapply(recepter, function(x){
      x = allSymbols[match(x,allSymbols)]
      x = na.omit(x)
      ans = bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
      return(ans[,2])
    })
    recepter_name = recepter_name[lengths(recepter)>0]
    recepter = recepter[lengths(recepter)>0]

    recepter2 = numeric(length(rank))
    #筛选含有相应受体的pathway
    pathway = c()
    n = c()
    s = 1
    results = kk@geneSets
    for (i in 1:length(results)) {

      a = sapply(recepter, function(x){
        ans = match(x,results[[i]])
        return(sum(!is.na(ans)) - length(x))

      })

      b = sum(test[,2]%in%results[[i]])
      if(any(a == 0) & b>minGenes){
        pathway[s] = i
        n[s] = b
        s = s+1
        c = which(a == 0)
        recepter2[c] = recepter2[c]+1
      }
    }
    if(length(pathway) == 0){
      next
    }
    use_recepter = recepter_name[recepter2>0 & lengths(recepter)>0]
    if(length(use_recepter) == 0){
      next
    }
    use_recepter = paste0(use_recepter,' ')
    R = c(R,use_recepter)
    recepter = recepter[recepter2>0 & lengths(recepter)>0]
    #构建nodes

    #受体nodes

    id = names(kk@geneSets)[pathway]
    nodes2 = as.data.frame(use_recepter)
    #下游通路nodes
    for (i in 1:length(pathway)) {
      b = keggGet(id[i])[[1]][['NAME']]
      b = str_split(b,'-')[[1]][1]
      nodes2 = rbind(nodes2,b)
      xiayou = rbind(xiayou,b)
    }
    colnames(nodes2) = 'name'
    nodes = as.data.frame(nodes[0:k+1,1])
    colnames(nodes) = 'name'
    go = as.data.frame(nodes[0:k,1])
    colnames(go) = 'name'
    shangyou = rbind(shangyou,go)
    nodes = rbind(nodes,nodes2)

    #links三列：出发，接受，连线粗细
    linkGenes = c()#连线到下游pathway的粗细
    s = 1
    results = kk@geneSets
    for (i in pathway) {
      a = test[,2]%in%results[[i]]
      linkGenes[s] = sum(a)
      s = s+1
    }

    #得到recepter到pathway的link
    link3 = c()
    s = 1
    nrec = length(use_recepter)
    for (i in pathway) {
      a = sapply(recepter, function(x){
        ans = match(x,results[[i]])
        return(sum(!is.na(ans)) - length(x))

      })
      b = which(a == 0)
      c = matrix(nrow = length(b),ncol = 4)
      c = as.data.frame(c)
      #赋值recepter序号
      c[,1] = b+k
      #赋值pathway序号
      c[,2] = k+nrec+s
      #赋值连线粗细
      c[,3] = linkGenes[s]/length(b)
      #赋值连接颜色
      c[,4] = color[b+1]
      link3 = rbind(link3,c)
      s = s+1
    }
    #得到ligand到recepter的link
    link2 = matrix(nrow = nrec,ncol = 4)
    link2 = as.data.frame(link2)
    #赋值ligand序号
    link2[,1] = k
    #赋值recepter序号
    link2[,2] = (k+1):(k+nrec)
    #赋值连线粗细
    s = 1
    for (i in 1:nrec) {
      a = which(link3[,1] == k+i)
      b = sum(link3[a,3])
      link2[i,3] = b
    }
    #赋值颜色，和link3中一致
    link2[,4] = color[2:(1+nrec)]
    #合并link
    colnames(link2) = c('source','target','value','color')
    colnames(link3) = c('source','target','value','color')
    links[0:k,'color'] = color[1]
    links = rbind(links[0:k,],link2,link3)
    z = sum(link2[,3])/k
    links[1:k,3] = z

    if(h == 1){
      links_co = links
      nodes_co = nodes
      links_co[,4] = color[p]
      p = p+1
      h = h+1
    }
    else{
      links_1 = links
      nodes_1 = nodes

      links_2 = links_co
      nodes_2 = nodes_co

      nodes_co = unique(rbind(nodes_1,nodes_2))
      co1 = match(nodes_1[,1],nodes_co[,1])
      co2 = match(nodes_2[,1],nodes_co[,1])

      #新ligand的links
      links_co1 = data.frame(matrix(nrow = nrow(links_1), ncol = 4))
      links_co1[,1] = co1[links_1[,1]+1]-1
      links_co1[,2] = co1[links_1[,2]+1]-1
      links_co1[,3] = links_1[,3]
      links_co1[,4] = color[p]
      p = p+1

      links_co2 = data.frame(matrix(nrow = nrow(links_2), ncol = 4))
      links_co2[,1] = co2[links_2[,1]+1]-1
      links_co2[,2] = co2[links_2[,2]+1]-1
      links_co2[,3] = links_2[,3]
      links_co2[,4] = links_2[,4]

      links_co = rbind(links_co1,links_co2)
      colnames(links_co) = c('source','target','value','color')
      colnames(nodes_co) = 'name'
    }

  }

  sankey_graph = list(links_co,nodes_co)
  return(sankey_graph)
}












































