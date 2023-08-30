###################################################################
## Load Mouse cell cycle genes
###################################################################
LoadMouseCellCycleGenes = function(){
  load("https://github.com/jginz/SCCodes/blob/main/mm2Hs.rda")
  g2m.features=as.character(na.exclude(mm2Hs[match(cc.genes.updated.2019$g2m.genes,mm2Hs$Human),"Mouse"]))
  s.features=as.character(na.exclude(mm2Hs[match(cc.genes.updated.2019$s.genes,mm2Hs$Human),"Mouse"]))
  return(list(g2m.features=g2m.features,s.features=s.features))
}

###################################################################
## Remove mitochondrial genes
###################################################################
RemoveMito = function(seu){
  remove_genes = rownames(seu)[grep("^mt-",rownames(seu))]
  seu = subset(seu,features=remove_genes,invert=T)
  return(seu)
}

###################################################################
## Remove RBC genes
###################################################################
RemoveRBC = function(seu){
  remove_genes = rownames(seu)[grep("^Hb[^(p)]",rownames(seu))]
  seu = subset(seu,features=remove_genes,invert=T)
  return(seu)
}

###################################################################
## Remove Ribosomal genes/proteins
###################################################################
RemoveRibo = function(seu){
  remove_genes = rownames(seu)[grep("^Rp[sl]",rownames(seu))]
  seu = subset(seu,features=remove_genes,invert=T)
  return(seu)
}

###################################################################
## Run doublet finder
##
## This script will need to be run on Seurat Object
## Please run PCA on the seurat object before running this script
###################################################################

RunDoubletFinder = function(seu){
  ## Find mini
  stdv = seu[["pca"]]@stdev
  sum.stdv = sum(seu[["pca"]]@stdev)
  percent.stdv = (stdv / sum.stdv) * 100
  cumulative = cumsum(percent.stdv)
  co1 = which(cumulative > 90 & percent.stdv < 5)[1]
  co2 = sort(which((percent.stdv[1:length(percent.stdv) - 1] - percent.stdv[2:length(percent.stdv)]) > 0.1),decreasing = T)[1]+1
  min.pc = min(co1, co2)
  
  ## Run Doublet finding algorithm
  sweep.list = paramSweep_v3(seu, PCs = 1:min.pc, num.cores = 2)
  sweep.stats = summarizeSweep(sweep.list)
  bcmvn = find.pK(sweep.stats)
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  ## Homotypic doublet proportion estimate
  annotations <- seu@meta.data$Clusterss
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(seu@meta.data)) 
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  # run DoubletFinder
  seu <- doubletFinder_v3(seu = seu,PCs = 1:min.pc,pK = optimal.pk,nExp = nExp.poi.adj)
  colnames(seu@meta.data)[grep("DF.classifications",names(head(seu)))]="doublet_finder"
  seu = subset(seu, doublet_finder == "Singlet")
  seu@meta.data = seu@meta.data[,colnames(seu@meta.data)[-grep("pANN_",colnames(seu@meta.data))]]
  seu = RunPCA(seu,verbose=F)
  return(seu)
}

###################################################################
## Find minimum PCs for downstream UMAP analysis
## Using ElbowPoint method
###################################################################

FindMinPCs = function(seu){
  feats = as.character(VariableFeatures(seu))
  tmp = subset(seu,features=feats)
  mat = as.matrix(tmp@assays$RNA@scale.data)
  p <- pca(mat, metadata = tmp@meta.data, removeVar = 0.1)
  elbow <- findElbowPoint(p$variance)
  return(as.numeric(elbow))
}

###################################################################
## UMAP themes (plotting function)
##
###################################################################
## UMAP theme
umap_theme <- theme(
  axis.line=element_blank(),
  axis.text.x=element_blank(),
  axis.text.y=element_blank(),
  axis.ticks=element_blank(),
  axis.title.x=element_blank(),
  axis.title.y=element_blank(),
  panel.background=element_blank(),
  panel.border=element_blank(),
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  plot.background=element_blank()
)

## Blank theme
theme_blank <- function(add_coord = TRUE, xlen_npc = 0.15, ylen_npc = 0.15, xlab = "", ylab = "", lab_size = 12, ...) {
  if (isTRUE(add_coord) && isTRUE(xlab != "")) {
    x_space <- lab_size + 2
  } else {
    x_space <- 0
  }
  if (isTRUE(add_coord) && isTRUE(ylab != "")) {
    y_space <- lab_size + 2
  } else {
    y_space <- 0
  }
  args1 <- list(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.background = element_blank(),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(0, 0, x_space, y_space, unit = "points"),
    complete = FALSE
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% formalArgs(theme)]
  out <- do.call(
    what = theme,
    args = args
  )
  if (isTRUE(add_coord)) {
    g <- grobTree(gList(
      linesGrob(x = unit(c(0, xlen_npc), "npc"), y = unit(c(0, 0), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
      textGrob(label = xlab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = 4 / 3, hjust = 0, gp = gpar(fontsize = lab_size)),
      linesGrob(x = unit(c(0, 0), "npc"), y = unit(c(0, ylen_npc), "npc"), arrow = arrow(length = unit(0.02, "npc")), gp = gpar(lwd = 2)),
      textGrob(label = ylab, x = unit(0, "npc"), y = unit(0, "npc"), vjust = -2 / 3, hjust = 0, rot = 90, gp = gpar(fontsize = lab_size))
    ))
    return(list(
      list(annotation_custom(g)),
      list(theme_scp() + out),
      list(coord_cartesian(clip = "off"))
    ))
  } else {
    return(list(
      list(theme_scp() + out)
    ))
  }
}

## SCP theme
theme_scp <- function(aspect.ratio = NULL, base_size = 12, ...) {
  text_size_scale <- base_size / 12
  args1 <- list(
    aspect.ratio = aspect.ratio,
    text = element_text(size = 12 * text_size_scale, color = "black"),
    plot.title = element_text(size = 14 * text_size_scale, colour = "black", vjust = 1),
    plot.subtitle = element_text(size = 13 * text_size_scale, hjust = 0, margin = margin(b = 3)),
    plot.background = element_rect(fill = "white", color = "white"),
    axis.line = element_blank(),
    axis.title = element_text(size = 13 * text_size_scale, colour = "black"),
    axis.text = element_text(size = 12 * text_size_scale, colour = "black"),
    strip.text = element_text(size = 12.5 * text_size_scale, colour = "black", hjust = 0.5, margin = margin(3, 3, 3, 3)),
    strip.background = element_rect(fill = "transparent", linetype = 0),
    strip.switch.pad.grid = unit(-1, "pt"),
    strip.switch.pad.wrap = unit(-1, "pt"),
    strip.placement = "outside",
    legend.title = element_text(size = 12 * text_size_scale, colour = "black", hjust = 0),
    legend.text = element_text(size = 11 * text_size_scale, colour = "black"),
    legend.key = element_rect(fill = "transparent", color = "transparent"),
    legend.background = element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    panel.border = element_rect(fill = "transparent", colour = "black", linewidth = 1),
    complete = TRUE
  )
  args2 <- as.list(match.call())[-1]
  call.envir <- parent.frame(1)
  args2 <- lapply(args2, function(arg) {
    if (is.symbol(arg)) {
      eval(arg, envir = call.envir)
    } else if (is.call(arg)) {
      eval(arg, envir = call.envir)
    } else {
      arg
    }
  })
  for (n in names(args2)) {
    args1[[n]] <- args2[[n]]
  }
  args <- args1[names(args1) %in% formalArgs(theme)]
  out <- do.call(
    what = theme,
    args = args
  )
  return(out)
}

###################################################################
## Run SCType
###################################################################
## Main function

#' @description GNU General Public License v3.0 
#' (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
#' @author Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

##' @custom_gene_list is the path to the gene list for SCType
##' @tissue_type is the unique character describing cell types
##' @plot plotting function for cell types

RunSCType <- function(seu,custom_gene_list=NULL,tissue_type="Immune system",plot=T){
  
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  
  if(is.null(custom_gene_list)){
    db_ = "https://github.com/jginz/SCCodes/blob/main/ScTypeDB_full.xlsx"
    gs_list = gene_sets_prepare(db_,tissue_type)
  }else{
    db_ = custom_gene_list
    gs_list = gene_sets_prepare(db_,tissue_type)
  }
  
  set.seed(1)
  es.max = sctype_score(as.matrix(GetAssayData(seu, slot = "data")), scaled = FALSE,
                        gs = gs_list$gs_positive)
  
  cL_results = do.call("rbind", lapply(unique(seu@meta.data$Clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seu@meta.data[seu@meta.data$Clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu@meta.data$Clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/6] = "Unknown"
  seu@meta.data$SCType = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seu@meta.data$SCType[seu@meta.data$Clusters == j] = as.character(cl_type$type[1])
  }
  if(plot){
    # prepare edges
    cL_results=cL_results[order(cL_results$cluster),]; 
    edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); 
    edges$cluster = paste0("cluster ", edges$cluster); 
    edges = edges[,c("cluster", "type")]; 
    colnames(edges) = c("from", "to"); 
    rownames(edges) <- NULL
    
    # prepare nodes
    nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; 
    nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); 
    nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; 
    nodes_lvl1$realname = nodes_lvl1$cluster; 
    nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
    ccolss = dittoColors()
    for (i in 1:length(unique(cL_results$cluster))){
      dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; 
      nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), 
                                                ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
    }
    nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
    files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")];
    files_db = unique(files_db); 
    nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
    nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; 
    nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
    
    mygraph <- graph_from_data_frame(edges, vertices=nodes)
    
    # Make the graph
    gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
      geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
      geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
      theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, 
                                        colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5))) + 
      geom_node_label(aes(filter=ord==1,label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), 
                      repel = !0, segment.linetype="dotted")
    gridExtra::grid.arrange(DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr)
    
  }
  return(seu)
}

###################################################################
## Run SCType
###################################################################
## Main function

#' @description GNU General Public License v3.0 
#' (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
#' @author Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

##' @custom_gene_list is the path to the gene list for SCType
##' @tissue_type is the unique character describing cell types
##' @plot plotting function for cell types

RunSCType2 = function(seu,custom_gene_list=NULL,tissue_type="Immune system",plot=T){
  
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  
  if(is.null(custom_gene_list)){
    db_ = "https://github.com/jginz/SCCodes/blob/main/ScTypeDB_full.xlsx"
    gs_list = gene_sets_prepare(db_,tissue_type)
  }else{
    db_ = custom_gene_list
    gs_list = gene_sets_prepare(db_,tissue_type)
  }
  
  set.seed(1)
  es.max = sctype_score(as.matrix(GetAssayData(seu, slot = "data")), scaled = FALSE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  cL_results = do.call("rbind", lapply(unique(seu@meta.data$Clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seu@meta.data[seu@meta.data$Clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu@meta.data$Clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  seu@meta.data$SCType2 = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seu@meta.data$SCType2[seu@meta.data$Clusters == j] = as.character(cl_type$type[1])
  }
  if(plot){
    # prepare edges
    cL_results=cL_results[order(cL_results$cluster),]; 
    edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); 
    edges$cluster = paste0("cluster ", edges$cluster); 
    edges = edges[,c("cluster", "type")]; 
    colnames(edges) = c("from", "to"); 
    rownames(edges) <- NULL
    
    # prepare nodes
    nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; 
    nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); 
    nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; 
    nodes_lvl1$realname = nodes_lvl1$cluster; 
    nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
    ccolss = dittoColors()
    for (i in 1:length(unique(cL_results$cluster))){
      dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; 
      nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), 
                                                ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
    }
    nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
    files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")];
    files_db = unique(files_db); 
    nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
    nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; 
    nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
    
    mygraph <- graph_from_data_frame(edges, vertices=nodes)
    
    # Make the graph
    gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
      geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
      geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
      theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, 
                                        colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5))) + 
      geom_node_label(aes(filter=ord==1,label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), 
                      repel = !0, segment.linetype="dotted")
    gridExtra::grid.arrange(DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr)
    
  }
  return(seu)
}

##################################################################
## Run SCType
###################################################################
## Main function

#' @description GNU General Public License v3.0 
#' (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
#' @author Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021

##' @custom_gene_list is the path to the gene list for SCType
##' @tissue_type is the unique character describing cell types
##' @plot plotting function for cell types

RunSCType3 = function(seu,custom_gene_list=NULL,tissue_type="Immune system",plot=T){
  
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  
  if(is.null(custom_gene_list)){
    db_ = "https://github.com/jginz/SCCodes/blob/main/ScTypeDB_full.xlsx"
    gs_list = gene_sets_prepare(db_,tissue_type)
  }else{
    db_ = custom_gene_list
    gs_list = gene_sets_prepare(db_,tissue_type)
  }
  
  set.seed(1)
  es.max = sctype_score(as.matrix(GetAssayData(seu, slot = "data")), scaled = FALSE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  cL_results = do.call("rbind", lapply(unique(seu@meta.data$Clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seu@meta.data[seu@meta.data$Clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu@meta.data$Clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  seu@meta.data$SCType3 = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seu@meta.data$SCType3[seu@meta.data$Clusters == j] = as.character(cl_type$type[1])
  }
  if(plot){
    # prepare edges
    cL_results=cL_results[order(cL_results$cluster),]; 
    edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); 
    edges$cluster = paste0("cluster ", edges$cluster); 
    edges = edges[,c("cluster", "type")]; 
    colnames(edges) = c("from", "to"); 
    rownames(edges) <- NULL
    
    # prepare nodes
    nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; 
    nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); 
    nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; 
    nodes_lvl1$realname = nodes_lvl1$cluster; 
    nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
    ccolss = dittoColors()
    for (i in 1:length(unique(cL_results$cluster))){
      dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; 
      nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), 
                                                ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
    }
    nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
    files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")];
    files_db = unique(files_db); 
    nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
    nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; 
    nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
    
    mygraph <- graph_from_data_frame(edges, vertices=nodes)
    
    # Make the graph
    gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
      geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + 
      geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
      theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, 
                                        colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5))) + 
      geom_node_label(aes(filter=ord==1,label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), 
                      repel = !0, segment.linetype="dotted")
    gridExtra::grid.arrange(DimPlot(seu, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr)
    
  }
  return(seu)
}

###################################################################
## Run Single Sample GSEA using UCell
## Source: 
## https://github.com/carmonalab/UCell
## https://www.sciencedirect.com/science/article/pii/S2001037021002816?via%3Dihub
###################################################################

## Geneset file is a custom .gmt file
## Default number of cores = 2
RunUCell = function(seu, geneset_file = "https://github.com/jginz/SCCodes/blob/main/acharya_genesets_mouse.gmt",ncores=2){
  gs = read.gmt(geneset_file)
  ES = enrichIt(obj = seu,gene.sets = gs, method = "UCell",groups = 1000, cores = ncores,min.size = 2)
  seu = AddMetaData(seu,ES)
  return(seu)
}

###################################################################
## Get pseudo-bulk RNAseq data based on a grouping variable
###################################################################

RunPseudoBulk = function(seu,grouping_variable="Clusters",return_as_Seurat=F){
  if(return_as_Seurat){
    out = AverageExpression(seu,group.by=grouping_variable,return.seurat = T)$RNA
  }else{
    out = data.frame(AverageExpression(seu,group.by=grouping_variable)$RNA,check.names=F)
  }
  return(out)
}
  


#################################################
## Percent of cells expressing gene of interest https://github.com/satijalab/seurat/issues/371#issuecomment-486384854
######################################################

PrctCellExpringGene <- function(object, genes, group.by = "all", assay = "RNA", datatype = "counts", threshold = 0){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object, assay = assay, datatype=datatype, threshold=threshold))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    results = lapply(list, PrctCellExpringGene, genes=genes, assay = assay, datatype=datatype, threshold=threshold)
    results %>% reduce(full_join, by="Markers") %>% select(any_of("Markers")) -> genelist
    results %>% reduce(full_join, by="Markers") %>% select(!any_of("Markers")) %>% "colnames<-"(names(results)) -> percentages
    combined <- cbind(genelist,percentages)
    return(combined)
  }
}

calc_helper <- function(object,genes,assay,datatype,threshold){
  counts = slot(object[[assay]],datatype)
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    round((sum(counts[genes,]>threshold)/ncells)*100,1)
  }else{return(NA)}
}


######## Error vector memory exhausted 
######### https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached


##################################################################
## Run Seurat pipeline on a single h5 file
#' @path path to .h5 file
###################################################################

RunSeuratPipeline = function(path=path_to_h5,for_integration=TRUE){
  dat = Read10X_h5(path_to_h5)
  seu = CreateSeuratObject(counts=dat)
  seu = PercentageFeatureSet(seu, "^mt-", col.name = "percent_mito")
  seu = PercentageFeatureSet(seu, "^Rp[sl]", col.name = "percent_ribo")
  seu = PercentageFeatureSet(seu, "^Hb[^(ap)]", col.name = "percent_hb")
  seu = subset(seu, subset = percent_mito<10 & percent_ribo > 5 & percent_hb < 1 )
  seu = subset(seu,subset = nCount_RNA>10 & nCount_RNA<20000 & nFeature_RNA > 10 & nFeature_RNA < 5000)
  seu = NormalizeData(seu,verbose=F)
  cc = LoadMouseCellCycleGenes()
  g2m_features = as.character(cc$g2m.features)
  s_features = as.character(cc$s.features)
  seu = CellCycleScoring(seu,g2m.features=g2m_features,s.features=s_features,verbose=F)
  seu = FindVariableFeatures(seu,selection.method="vst",nfeatures = 2000,verbose=F)
  seu = ScaleData(seu,vars.to.regress=c("nCount_RNA","S.Score","G2M.Score","percent_mito","percent_hb","percent_ribo"),verbose=F)
  seu = RunPCA(seu,verbose=F)
  seu = RunDoubletFinder(seu)
  if(for_integration==FALSE){
    seu = min_pcs = ifelse(max(FindMinPCs(seu),10)>15,10,10)
    seu = RunUMAP(seu, dims = 1:min_pcs,umap.method="uwot",n.neighbors = 50,min.dist=0.01,metric="cosine",n.components=2,verbose=F)
    seu = FindNeighbors(object = seu, dims = 1:min_pcs,verbose=F)
    seu = FindClusters(object = seu, resolution = c(0.4,0.6,0.8),verbose=F,save.SNN = T)
    Final_Resolution = "RNA_snn_res.0.6"  ## Identify your final cluster resolution
    Idents(seu)<-Final_Resolution
    seu@meta.data$Clusters = factor(as.character(as.numeric(Idents(seu))),levels=1:length(Idents(seu)))
    Idents(seu)<-'Clusters'
    seu = RunSCType(seu)
  }
  return(seu)
}

###################################################################
## Integrate data using Seurat pipeline
###################################################################

RunSeuratIntegration = function(h5_file_list){
  seu_list = lapply(h5_file_list, RunSeuratPipeline)
  seu_features = SelectIntegrationFeatures(object.list = seu_list, nfeatures = 2000)
  anchors = FindIntegrationAnchors(object.list = seu_list, scale=F, normalization.method = "LogNormalize", anchor.features = seu_features,dims = 1:30)
  dat_int = IntegrateData(anchorset = anchors, dims = 1:30)
  dat_int = ScaleData(dat_int, verbose = FALSE)
  dat_int = RunPCA(object = dat_int, verbose = FALSE)
  dat_int = RunUMAP(object = dat_int, dims = 1:10, umap.method="uwot",n.neighbors = 10L,min.dist=0.01,metric="cosine",verbose=F)
  dat_int = FindNeighbors(dat_int, reduction = "pca", dims = 1:10)
  dat_int = FindClusters(dat_int,resolution=c(0.4,0.6,0.8),random.seed=123)
  return(dat_int)
}
