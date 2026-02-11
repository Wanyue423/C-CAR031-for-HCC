library(Seurat)
library(tidyverse)

source("code base/spatial analysis/spatial analysis.R")
source("code base/color_pallete.R")

out_dir_v5 <- "20251022 epi nn analysis"
dir.create(out_dir_v5)

all_data <- readRDS("report_all_samples/output/data/adata_celltype.rds")

# define function ------
#### flatDoubleList
flatDoubleList <- function(x, sep = "---"){
  stopifnot(is.vector(x))
  
  res_ls <- list()
  
  # layer 1
  name_1 <- names(x)
  id_1   <- seq_along(x)
  if(is.null(name_1)) name_1 <- id_1
  
  for(i in id_1){
    y <- x[[i]]
    stopifnot(is.vector(y))
    
    # layer 2
    name_2 <- names(y)
    id_2   <- seq_along(y)
    if(is.null(name_2)) name_2 <- id_2
    
    for(j in id_2){
      res  <- y[[j]]
      name <- paste0(name_1[i], sep, name_2[j])
      res_ls[[name]] <- res
    }
  }
  res_ls
}


# cn analysis ------------
cells_position_df <-
  Embeddings(all_data[['spatial']]) %>%
  as.data.frame() %>%
  rownames_to_column('cell_id') %>%
  mutate(
    x_centroid = SPATIAL_1,
    y_centroid = SPATIAL_2,
    cell_id2 = cell_id
  ) %>%
  column_to_rownames("cell_id2")

all_data_meta_for_cn_analysis <- 
  # all_data[[c("project", "group", "annotation_full2")]] %>% 
  # dplyr::rename(annotation_full = annotation_full2 ) %>% 
  
  all_data[[]] %>%
  mutate(celltype = as.character(celltype)) %>% 
  mutate(project = sample) %>%
  mutate(annotation = ifelse(celltype %in% t_cells, "T", celltype)) %>%
  
  dplyr::filter(
    ! annotation %in% c("low-quality cell", "Hepatocyte")
  ) %>% 
  
  dplyr::mutate(
    annotation_full = annotation,
    cell_id_raw = cell_id
  ) %>% 
  
  dplyr::select(-cell_id) %>% 
  
  rownames_to_column("cell_id") %>%
  left_join(cells_position_df[c('cell_id', 'x_centroid', 'y_centroid')], by = "cell_id") %>%
  rename(X = x_centroid, Y = y_centroid) %>%
  mutate(
    annotation_full = factor(annotation_full), # import!
    obj_id = cell_id # import!
  )

# res_cn <- CN_compute_CN2(
#   all_data_meta_for_cn_analysis,
#   
#   query =  c("Tumor cell", "Proliferating tumor cell"),
#   
#   x = "X",
#   y = "Y",
#   group.by       = "project",
#   classification = "annotation_full",
#   return.vars    =  "freq", # "count",
#   k.nn           = 100, # 80um radius in papers ( ~100 cells )
#   k.cluster      = 5
# )

# res_cn <- readRDS(file.path("CN_analysis/res_cn.rds"))

# cluster by cell prop table --------------
all_data_cn <-
  CreateSeuratObject(
    CreateAssay5Object(
      counts = res_cn$cn_nn_cellfreq,
      data = res_cn$cn_nn_cellfreq
    ),
    meta.data = all_data_meta_for_cn_analysis[colnames(res_cn$cn_nn_cellfreq), ]
  )

all_data_cn <-
  all_data_cn %>% 
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(all_data_cn)

all_data_cn <- 
  all_data_cn %>%
  RunUMAP(dims = 1:5) %>%
  FindNeighbors(dims = 1:5, k.param = 20) %>% 
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5))

# saveRDS(all_data_cn, file.path(out_dir_v5, 'all_data_cn.rds'))


DimPlot(all_data_cn, group.by = "RNA_snn_res.0.1", label = TRUE) + NoLegend() 

pdf(
  file = file.path(out_dir_v5, "cn cluster - tumor cluster prop of group - add pre.pdf"),
  height = 6,
  width  = 4
)
all_data_cn[[]] %>% 
  mutate(tumor = str_c('Tumor C', as.character(RNA_snn_res.0.1), sep = "")) %>% 
  
  # dplyr::filter(group != "pre") %>%
  
  group_by(group, tumor) %>%
  summarise(count = n()) %>% 
  mutate(freq = count / sum(count))  %>%
  
  ggplot(aes(x = group, y = freq, fill = tumor)) +
  geom_col(width = 0.55) +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(
    x = NULL,
    y = NULL
  )
dev.off()

# count of nn cells ----------------
res_nn_ls <-
  all_data_meta_for_cn_analysis %>%
  group_by(sample) %>%
  nest() %>%
  deframe() %>%
  
  imap(function(dat, sample_name){
    query_dat <- dat[dat$annotation_full %in% c("Tumor cell", "Proliferating tumor cell"), ]
    nn <- 
      RANN::nn2(
        dat[c('X', 'Y')],
        query = query_dat[c('X', 'Y')],
        k = 100
      )
    nn$dat <-  dat
    nn$dat_idx <- dat$obj_id
    nn$query_idx <- query_dat$obj_id
    nn
  })

# saveRDS(res_nn_ls, file.path(out_dir_v5, 'res_nn_ls.rds'))

res_nn_name_ls <-
  res_nn_ls %>%
  map(function(x){ # browser()
    x$nn.idx %>% dim()
    
    x$nn.idx[1:20, 1:5]
    all_idx <- c(x$nn.idx) + 1
    all_celltype <- c("Others", as.character(x$dat$annotation_full))
    
    all_idx2celltype <- all_celltype[all_idx]
    all_idx2celltype_mat <-  
      matrix(
        all_idx2celltype,
        ncol = 100,
        byrow = FALSE
      )
    rownames(all_idx2celltype_mat) <- x$query_idx
    all_idx2celltype_mat
  })

res_nn_name_df <- 
  do.call(rbind, res_nn_name_ls) %>%
  as.data.frame() %>%
  rownames_to_column('obj_id') %>%
  left_join(
    all_data_meta_for_cn_analysis [c('celltype', 'project', 'group', 'annotation', 'annotation_full', 'obj_id')],
    by = "obj_id"
  ) %>%
  left_join(
    all_data_cn[[c('RNA_snn_res.0.1', 'obj_id')]],
    by = "obj_id"
  ) %>%
  mutate(tumor = str_c('Tumor C' , as.character(RNA_snn_res.0.1), sep = ""))

res_nn_name_df_full <-
  res_nn_name_df %>%
  pivot_longer(
    V1:V100,
    names_to = "NN_index",
    values_to = "NN_cells"
  ) 

# GPC3 cells  -------------- 
gpc_de_ls <-
  list.files('.', pattern = "GPC3") %>%
  setNames(., str_remove(., "\\..*$")) %>%
  map(readxl::read_excel)

gpc_de_gene_ls <-
  gpc_de_ls %>%
  map(function(x){
    de <- x %>%
      dplyr::filter(
        logFC > 1  & adj.P.Val < 0.05
      )
    
    # return
    list(
      gene = de$symbol,
      gene_id = de$ENTREZID
    )
  })
  
all_data_epi %>%
  FetchData(vars = gpc_de_gene_ls$GPC3_GROUP$gene) %>%
  cor %>%
  corrplot::corrplot(
    order = "hclust",
    col = corrplot::COL2('RdBu', 200) %>% rev()
  ) 
  
all_data_epi %>%
  FetchData(vars = gpc_de_gene_ls$GPC3_GROUP$gene) %>%
  colMeans() %>%
  enframe() %>%
  ggplot(aes(x = reorder(name, value), y = value)) +
  geom_col() +
  RotatedAxis() +
  geom_hline(yintercept = 0.1)



gpc_de_ls$GPC3_GROUP %>% 
  mutate(group = case_when(
    logFC > 1  & adj.P.Val < 0.05 ~ 'up',
    logFC < -1 & adj.P.Val < 0.05 ~ 'down',
    TRUE ~ 'non'
  )) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
  geom_point() +
  geom_text(
    data = ~ dplyr::filter(., symbol == 'GPC3'),
    aes(label = symbol)
  )

gpc_de_ls$GPC3_ALL %>% 
  mutate(group = case_when(
    logFC > 1  & adj.P.Val < 0.05 ~ 'up',
    logFC < -1 & adj.P.Val < 0.05 ~ 'down',
    TRUE ~ 'non'
  )) %>% 
  ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = group)) +
  geom_point() +
  geom_text(
    data = ~ dplyr::filter(., symbol == 'GPC3'),
    aes(label = symbol)
  )

# all_data_epi
all_data_epi <- all_data[, all_data$cell_id_sample %in% all_data_cn$cell_id_sample] 
all_data_epi <- RenameCells(all_data_epi, new.names = all_data_epi$cell_id_sample)


all_data_epi@meta.data <-
  all_data_epi@meta.data %>%
  rownames_to_column('id') %>%
  left_join(
    all_data_cn[[c(
      'cell_id_sample',
      'X', 'Y', 
      'RNA_snn_res.0.1', 
      'RNA_snn_res.0.2', 
      'RNA_snn_res.0.3',
      'RNA_snn_res.0.4', 
      'RNA_snn_res.0.5'
    )]],
    by = 'cell_id_sample'
  ) %>%
  column_to_rownames('id')
  
all_data_epi$annotation_epi <- str_c('Tumor C', all_data_epi$RNA_snn_res.0.1, sep = "")


## AddMouduleScore - gpc3 score -------
all_data_epi_gpc3 <- 
  all_data_epi %>%
  AddModuleScore(
    map(gpc_de_gene_ls, 'gene')
  )

# > gpc_de_gene_ls %>% names
# [1] "GPC3_ALL"   "GPC3_GROUP"


VlnPlot(
  all_data_epi_gpc3,
  features = 'Cluster2',
  pt.size = 0,
  group.by = 'RNA_snn_res.0.1'
)


# export
pdf(
  file = file.path(out_dir_v5, "module score - vlnplot.pdf"),
  height = 6,
  width = 8
)
  VlnPlot(
    all_data_epi_gpc3,
    features = 'Cluster1',
    pt.size = 0,
    group.by = 'RNA_snn_res.0.1'
  ) + ggtitle('GPC3_ALL')
  VlnPlot(
    all_data_epi_gpc3,
    features = 'Cluster2',
    pt.size = 0,
    group.by = 'RNA_snn_res.0.1'
  ) + ggtitle('GPC3_GROUP')
dev.off()


# dimplot
all_data_cn <- 
  AddMetaData(
    all_data_cn,
    all_data_epi_gpc3[[c("Cluster1", "Cluster2")]]
  )

all_data_cn2 <- all_data_cn

# all_data_cn2$Cluster1 <- NULL
# all_data_cn2$Cluster2 <- NULL

all_data_cn2@meta.data <-
  all_data_cn2@meta.data %>%
  rownames_to_column('id') %>%
  left_join(
    all_data_epi_gpc3[[c('cell_id_sample',"Cluster1", "Cluster2")]],
    by = 'cell_id_sample'
  ) %>%
  column_to_rownames('id')
 

# export
pdf(
  file = file.path(out_dir_v5, "module score - feature plot.pdf"),
  height = 6,
  width = 8
)
  all_data_cn2 %>%
    FeaturePlot(
      features = "Cluster1"
    ) + ggtitle('GPC3_ALL')
  all_data_cn2 %>%
    FeaturePlot(
      features = "Cluster2"
    ) + ggtitle('GPC3_GROUP')
dev.off()



##  AUCell ----------
# library(AUCell)

exprMatrix <- all_data_epi@assays$RNA@data
cell_ranking <- AUCell::AUCell_buildRankings(exprMatrix)

cells_AUC <- AUCell::AUCell_calcAUC(map(gpc_de_gene_ls, 'gene'), cell_ranking)
set.seed(1234)

pdf(
  file.path(out_dir_v5, 'AUCell - GPC3 signature.pdf'),
  height = 6,
  width = 8
)
cells_assignment <- AUCell::AUCell_exploreThresholds(cells_AUC, plotHist = T, assignCells = T, smallestPopPercent = 0.1)
dev.off()


auc_thr1 <- cells_assignment$GPC3_GROUP$aucThr$thresholds['R_k3', 'threshold']
auc_thr1

all_data_epi$AUC_GPC3 <- as.numeric(AUCell::getAUC(cells_AUC_TLS))
all_data_epi$AUC_GPC3_Group <- ifelse(all_data_epi$AUC_GPC3 > auc_thr1, yes = "GPC3", no = "Others")

p_gpc3_aucell <-
  all_data_epi[[]] %>%
  group_by(
    annotation_epi, 
    AUC_GPC3_Group
  ) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_epi, y = freq, fill = AUC_GPC3_Group )) +
  geom_col() +
  theme_bw() +
  ggsci::scale_fill_aaas()

ggsave(
  p_gpc3_aucell,
  filename = file.path(out_dir_v5, "AUCell - tumor cluster props.pdf"),
  height = 6,
  width = 8
)

# export VlnPlot 
pdf(
  file = file.path(out_dir_v5, "AUCell - vlnplot.pdf"),
  height = 6,
  width = 8
)
VlnPlot(
  all_data_epi,
  features = 'AUC_GPC3',
  pt.size = 0,
  group.by = 'annotation_epi'
) 
dev.off()

# de analysis and enrichment analysis --------
# all_data_epi$annotation_epi <- str_c('Tumor C', all_data_epi$RNA_snn_res.0.1, sep = "")
Idents(all_data_epi) <- 'annotation_epi'

all_marker_epi <-
  FindAllMarkers(
    all_data_epi,
    logfc.threshold = 0,
    return.thresh = 1.1
  )

  
# expot enrichment table
write.csv(
  all_marker_epi,
  file = file.path(out_dir_v5, str_glue("de analysis - de table.csv")),
  row.names = TRUE
)
  


# saveRDS(all_data_epi, 'all_data_epi.rds')


# enrichment analysis -----------
all_gene_db <- 
  clusterProfiler::bitr(
    rownames(all_data_epi),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  )

all_geneid_epi_ls <-
  list(epi = all_marker_epi) %>%
  map(function(dat){
    gene <-
      dat %>%
      dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05, pct.1 > 0.1) %>%
      group_by(cluster) %>%
      nest %>%
      deframe() %>%
      map(~{
        clusterProfiler::bitr(
          .x$gene,
          fromType = "SYMBOL",
          toType   = "ENTREZID",
          OrgDb    = "org.Hs.eg.db"
        ) %>%
          .$ENTREZID
      })
    
    gene  
  })

all_geneid_epi_ls <- 
  all_geneid_epi_ls$epi

all_enrich_epi_ls <- 
  all_geneid_epi_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      universe = all_gene_db$ENTREZID,
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment 
all_enrich_epi_ls <-
  all_enrich_epi_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    all_genes <- unlist(xx)
    gene_db <- 
      clusterProfiler::bitr(
        all_genes,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = "org.Hs.eg.db"
      ) 
    
    xx2 <-
      xx %>%
      map_chr(function(g){
        str_c(gene_db$SYMBOL[match(g, gene_db$ENTREZID)], collapse = ",")
      })
    x@result$geneSymbol <- xx2
    x
  })


# enrichment plot
p_all_enrich_epi_ls <-
  all_enrich_epi_ls %>%
  imap(function(enrich, name){ # browser()
    enrichplot::dotplot(
      enrich,
      x = "GeneRatio",
      color = "p.adjust",
      showCategory = 10,
      label_format = 40
    ) + 
      # scale_fill_gradient2(low = "gray", high = "red") +
      ggtitle(name)
  })


# expot enrichment table
all_enrich_epi_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(out_dir_v5, str_glue("de enrichment analysis - enriment table.csv")),
      row.names = FALSE
    )
  })

# expot enrichment plot
pdf(
  file.path(out_dir_v5, "de enrichment analysis dotplot.pdf"), 
  height = 6.5, 
  width  = 8
)
p_all_enrich_epi_ls %>%  walk(print)
dev.off()


#  msigdb hallmark -----
h_geneset <- msigdbr::msigdbr(category = "H")

h_geneset$gs_name %>%
  unique() %>% 
  str_to_lower() %>%
  str_remove("^hallmark_")

h_geneset_ls <-
  h_geneset[c('gs_name', 'gene_symbol')] %>%
  with(split(gene_symbol, gs_name)) %>%
  {set_names(
    .,
    names(.) %>% 
      str_to_lower() %>%
      str_remove("^hallmark_")
  )}

all_data_epi_score <-
  all_data_epi %>%
  AddModuleScore(
    features = h_geneset_ls,
    # name = names(h_geneset_ls)
  )

all_data_epi_score %>% head

# msigdb anno
h_score_anno <- 
  all_data_epi_score[[]] %>%
  dplyr::select(annotation_epi, Cluster1:Cluster50) %>%
  group_by(annotation_epi) %>%
  summarise(across(everything(), mean))

p_hallmark_1 <-
  pheatmap::pheatmap(
    h_score_anno %>% column_to_rownames("annotation_epi") %>% t,
    scale = 'row',
    labels_row = names(h_geneset_ls),
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
    border_color = NA,
    silent = TRUE
  )

# export 
pdf(
  file.path(out_dir_v5, "Hallmark pathway score - tumor cluster.pdf"), 
  height = 12, 
  width  = 5
)
p_hallmark_1$gtable %>% cowplot::plot_grid() %>% print()
dev.off()


# cluster  of epithelial cells -----------
epi_out_dir <- "re-cluster of epithelial"
dir.create(epi_out_dir)


all_data_raw <- readRDS("report_all_samples/output/data/adata_raw.rds")
all_data_raw <- RenameCells(all_data_raw, new.names = all_data_raw$cell_id_sample)

all_data_epi2 <- all_data_raw[, colnames(all_data_epi)]
all(colnames(all_data_epi) %in% colnames(all_data_epi2))

# only used PD_SD
# all_data_epi2 <- all_data_epi2[, all_data_epi2$group == "PD_SD"]

all_data_epi2 <- 
  all_data_epi2 %>%
  AddMetaData(
    all_data_epi[[c('AUC_GPC3', 'AUC_GPC3_Group')]]
  )

all_data_epi2$project <- all_data_epi2$sample
all_data_epi2$nCount_RNA <- all_data_epi2$total_counts
all_data_epi2$cell_id_ori <- all_data_epi2$cell_id
all_data_epi2$cell_id <- NULL



DefaultAssay(all_data_epi2) <- "RNA"


all_data_epi2 <-
  all_data_epi2 %>%
  NormalizeData(scale.factor = median(.$nCount_RNA)) %>%
  FindVariableFeatures(nfeatures = 2000) 

VariableFeatures(all_data_epi2) <- gpc_de_gene_ls$GPC3_GROUP$gene
# VariableFeaturePlot(all_data_epi2)
all_data_epi2 <- all_data_epi2 %>% ScaleData(vars.to.regress = 'nCount_RNA')



## pca -----
all_data_epi2 <- all_data_epi2 %>% RunPCA()

pdf(file.path(epi_out_dir, 'pca plot.pdf'), height = 6, width = 8)
ElbowPlot(all_data_epi2, ndims = 50)
dev.off()

set.seed(1234)
suset_cells_epi <- 
  all_data_epi2[[]] %>%
  rownames_to_column("cell_id") %>%
  group_by(project) %>%
  slice_sample(prop = 50000/ncol(all_data_epi2) ) %>%
  dplyr::pull("cell_id")


pdf(file.path(epi_out_dir, 'dimheatmap 1-5.pdf'), height = 6, width = 10)
DimHeatmap(all_data_epi2[, suset_cells_epi], dims = 1:5)
dev.off()
pdf(file.path(epi_out_dir, 'dimheatmap 6-10.pdf'), height = 6, width = 10)
DimHeatmap(all_data_epi2[, suset_cells_epi], dims = 6:10)
dev.off()
pdf(file.path(epi_out_dir, 'dimheatmap 11-15.pdf'), height = 6, width = 10)
DimHeatmap(all_data_epi2[, suset_cells_epi], dims = 11:15)
dev.off()
pdf(file.path(epi_out_dir, 'dimheatmap 16-20.pdf'), height = 6, width = 10)
DimHeatmap(all_data_epi2[, suset_cells_epi], dims = 16:20)
dev.off()


## dimension reduction --------
epi_dims <- 15 #  10 # 20

# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb

all_data_epi2 <- RunUMAP(all_data_epi2, dims = 1:epi_dims)


pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} project.pdf")), height = 6, width = 8)
DimPlot(all_data_epi2, group.by = "project", label = TRUE) + NoLegend()
dev.off()


all_data_epi2 <- all_data_epi2 %>% 
  harmony::RunHarmony(
    dims.use = 1:epi_dims,
    group.by.vars = "project",
    reduction.save = "harmony"
    # theta = 3, # default 2, bigger value -- more diversity per cluster
    # lambda = 0.5, # default, smaller value -- more integrated
  )

all_data_epi2 <- RunUMAP(
  all_data_epi2, 
  dims = 1:epi_dims, 
  reduction = "harmony", 
  reduction.name = "harmonyumap", 
  reduction.key  = "harmonyUMAP_",
  return.model = TRUE
)


pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} project - harmony.pdf")), height = 6, width = 8)
DimPlot(all_data_epi2, group.by = "project", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()


reduction.use <- 'harmony' # 'harmony' 'pca'
all_data_epi2 <-
  FindNeighbors(
    all_data_epi2,
    reduction = reduction.use,
    dims = 1:epi_dims,
    k.param = 10
  )

all_data_epi2 <-
  FindClusters(
    all_data_epi2,
    resolution = c(0.05, 0.1, 0.2) #, 0.4, 0.6, 0.8)
  )

# dimplot
pdf(file.path(epi_out_dir, str_glue("Dimplot dim{epi_dims} cluster - {reduction.use}.pdf")), height = 6, width = 8)
DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "umap") + NoLegend()
# DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "umap") + NoLegend()

DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.05", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.1", label = TRUE, reduction = "harmonyumap") + NoLegend()
DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.2", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.4", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.6", label = TRUE, reduction = "harmonyumap") + NoLegend()
# DimPlot(all_data_epi2, group.by = "RNA_snn_res.0.8", label = TRUE, reduction = "harmonyumap") + NoLegend()
dev.off()

saveRDS(all_data_epi2, file = file.path(epi_out_dir, "all_data_epi2.rds"))


# dimplot 
DimPlot(
  all_data_epi2, 
  group.by = "AUC_GPC3_Group",
  reduction = "harmonyumap"
)

all_data_epi2[[]] %>%
  group_by(group, AUC_GPC3_Group, RNA_snn_res.0.2) %>%
  summarise(count = n()) %>%
  group_by(group, AUC_GPC3_Group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = AUC_GPC3_Group, y = freq, fill = RNA_snn_res.0.2)) +
  facet_wrap(~ group) +
  geom_col()
  
all_data_epi2[[]] %>%
  group_by(group, RNA_snn_res.0.2, AUC_GPC3_Group) %>%
  summarise(count = n()) %>%
  group_by(group, RNA_snn_res.0.2) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = RNA_snn_res.0.2, y = freq, fill = AUC_GPC3_Group)) +
  facet_wrap(~ group, scales = "free_x") +
  geom_col()


all_data_epi2 %>%
  VlnPlot(
    group.by = "RNA_snn_res.0.2",
    features = "AUC_GPC3",
    pt.size = 0
  )



# 
# all_data_epi2 %>%
#   DimPlot(
#     
#   )


# GPC3 negative cell -------------

# total tumor props, PD > PR
# PD中GPC3- prop > GPC3 + 


# DSP GP3 gene 


##  gpc3 cell props   -------------
p_col <-
  all_data_epi2[[]] %>%  
  group_by(group, AUC_GPC3_Group) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq, fill = AUC_GPC3_Group )) +
  geom_col()

# export 
pdf(
  file.path(out_dir_v5, "GPC3 cells - cell prop of group.pdf"), 
  height = 8, 
  width  = 5
)
p_col %>% cowplot::plot_grid() %>% print()
dev.off()


## de analysis   -------------
Idents(all_data_epi2) <- "AUC_GPC3_Group" 
all_marker_gpc3_netg <-
  FindMarkers(
    all_data_epi2,
    ident.1 = "PR",
    ident.2 = "PD_SD",
    group.by = "group",
    subset.ident = "Others"
  )

de_df_gpc3 <-
  all_marker_gpc3_netg %>%
  mutate(
    group = case_when(
      avg_log2FC > 0.25 & p_val_adj < 0.05 ~ 'Pos',
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ 'Neg',
      TRUE ~ "Non"
    )
  ) 

write.csv(
  de_df_gpc3, 
  file = file.path(out_dir_v5, "GPC3 cells - de analysis table of PR vs PD_SD in others.csv"),
  row.names = TRUE
)


all_marker_gpc3_pos <-
  FindMarkers(
    all_data_epi2,
    ident.1 = "PR",
    ident.2 = "PD_SD",
    group.by = "group",
    subset.ident = "GPC3"
  )
all_marker_gpc3_pos %>%
  mutate(
    group = case_when(
      avg_log2FC > 0.25 & p_val_adj < 0.05 ~ 'Pos',
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ 'Neg',
      TRUE ~ "Non"
    )
  ) %>%
  write.csv(
    file = file.path(out_dir_v5, "GPC3 cells - de analysis table of PR vs PD_SD in GPC3.csv"),
    row.names = TRUE
  )
# rownames(de_df_gpc3)[de_df_gpc3$group == "Neg"]

p_volcano <-
  de_df_gpc3 %>%
  ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
  geom_point() +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = 2) 

  
  scale_color_manual(values = c(
    Neg = "blue",
    Pos = "red",
    Non = "gray"
  )) +
  ggtitle("PR vs PD_SD")


pdf(
  file.path(out_dir_v5, "GPC3 cells - volcano plot of GPC3 epi cells vs others epi cells.pdf"), 
  height = 5, 
  width  = 8
)
p_volcano %>% print()
dev.off()



## exporess of de gene  -------------

pdf(
  file.path(out_dir_v5, "GPC3 cells - dotplot of PD_SD up-regulation gene of others epi cells.pdf"), 
  height = 8, 
  width  = 5
)
DotPlot(
  all_data_epi2[, all_data_epi2$AUC_GPC3_Group == "Others"],
  features = rownames(de_df_gpc3)[de_df_gpc3$group == "Neg"],
  group.by = "group"
) +
  coord_flip()
dev.off()




# h_score_anno2 <- 
#   all_data_epi_score[[]] %>%
#   dplyr::select(group, AUC_GPC3_Group, Cluster1:Cluster50) %>%
#   group_by(group, AUC_GPC3_Group) %>%
#   summarise(across(everything(), mean)) %>%
#   mutate(id = str_c(group, AUC_GPC3_Group, sep = "_")) %>%
#   ungroup() %>%
#   dplyr::select(-c(group, AUC_GPC3_Group))
# 
# 
# p_hallmark_2 <-
#   pheatmap::pheatmap(
#     h_score_anno2 %>% column_to_rownames("id") %>% t,
#     scale = 'row',
#     labels_row = names(h_geneset_ls),
#     color = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
#     border_color = NA,
#     silent = TRUE
#   )

# cowplot::plot_grid(p_hallmark_2$gtable)

# 
# h_geneset$gs_name %>% unique
# 
# HALLMARK_PI3K_AKT_MTOR_SIGNALING
# 



## dotplot of PI3K-AKT gene  -------------
pdf(
  file.path(out_dir_v5, "GPC3 cells - dotplot of HALLMARK_PI3K_AKT_MTOR_SIGNALING gene.pdf"), 
  height = 12, 
  width  = 5
)
DotPlot(
  all_data_epi2[, all_data_epi2$AUC_GPC3_Group == "Others"],
  features = h_geneset$gene_symbol[h_geneset$gs_name == "HALLMARK_PI3K_AKT_MTOR_SIGNALING"],
  group.by = "group"
) +
  coord_flip()
dev.off()


## TGFb expresion ----------
# VlnPlot(
#   all_data_epi2[, all_data_epi2$AUC_GPC3_Group == "Others"],
#   features = c("TGFB1", "TGFB1I1", "TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3"),
#   group.by = "group",
# )

p_voln_tgfb_grp <-
  all_data_epi2 %>%
  FetchData(
    vars = c("TGFB1", "TGFB1I1", "TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3") %>% c("group", "AUC_GPC3_Group")
  ) %>%
  dplyr::filter(AUC_GPC3_Group == "Others") %>%
  pivot_longer(
    -c("group", "AUC_GPC3_Group"),
    names_to = "gene",
    values_to = "value"
  ) %>%
  dplyr::filter(value != 0) %>%
  
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_violin(trim = F) +
  ggpubr::geom_pwc(tip.length = 0) +
  facet_wrap(~ gene ) +
  theme_bw()

# 20251215 export 
p_voln_tgfb_grp$data %>% 
  write.csv(
    file = file.path(
      out_dir_v5, 
      "GPC3 cells - vonplot of TGFb gene in Others epi cells - group - remvoe zero values.csv"
    )
  )
### 




p_voln_tgfb_gpc3 <-
  all_data_epi2 %>%
  FetchData(
    vars = c("TGFB1", "TGFB1I1", "TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3") %>% c("group", "AUC_GPC3_Group")
  ) %>%
  # dplyr::filter(AUC_GPC3_Group == "Others") %>%
  pivot_longer(
    -c("group", "AUC_GPC3_Group"),
    names_to = "gene",
    values_to = "value"
  ) %>%
  dplyr::filter(value != 0) %>%
  
  ggplot(aes(x = AUC_GPC3_Group, y = value, fill = AUC_GPC3_Group)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.15, outliers = F) +
  ggpubr::geom_pwc(tip.length = 0) +
  facet_grid(group ~ gene) +
  theme_bw()
p_voln_tgfb_gpc3

# DotPlot(
#   all_data_epi2[, all_data_epi2$AUC_GPC3_Group == "Others"],
#   features = c("TGFB1", "TGFB1I1", "TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3"),
#   group.by = "group"
# ) +
#   coord_flip()


# export
pdf(
  file.path(out_dir_v5, "GPC3 cells - vonplot of TGFb gene in Others epi cells - group - remvoe zero values.pdf"), 
  height = 10, 
  width  = 12
)
p_voln_tgfb_grp %>% print()
dev.off()

pdf(
  file.path(out_dir_v5, "GPC3 cells - vonplot of TGFb gene in Others epi cells - group and GPC3 group - remvoe zero values.pdf"), 
  height = 10, 
  width  = 20
)
p_voln_tgfb_gpc3 %>% print()
dev.off()


# 20251215 
p_voln_tgfb_grp_gpc3 <-
  all_data_epi2 %>%
  FetchData(
    vars = c("TGFB1", "TGFB1I1", "TGFB2","TGFB3","TGFBR1","TGFBR2","TGFBR3") %>% c("group", "AUC_GPC3_Group")
  ) %>%
  dplyr::filter(AUC_GPC3_Group == "GPC3") %>%
  pivot_longer(
    -c("group", "AUC_GPC3_Group"),
    names_to = "gene",
    values_to = "value"
  ) %>%
  dplyr::filter(value != 0) %>%
  
  ggplot(aes(x = group, y = value, fill = group)) +
  geom_violin(trim = F) +
  ggpubr::geom_pwc(tip.length = 0) +
  facet_wrap(~ gene ) +
  theme_bw()


# export
pdf(
  file.path(out_dir_v5, "GPC3 cells - vonplot of TGFb gene in GPC3 epi cells - group - remvoe zero values.pdf"), 
  height = 10, 
  width  = 12
)
p_voln_tgfb_grp_gpc3 %>% print()
dev.off()


p_voln_tgfb_grp_gpc3$data %>% 
  write.csv(
    file = file.path(
      out_dir_v5, 
      "GPC3 cells - vonplot of TGFb gene in GPC3 epi cells - group - remvoe zero values.csv"
    )
  )
### 


## GSEA: KEGG and GO -----------
all_gene_df <- 
  clusterProfiler::bitr(
    rownames(all_data_epi),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ) %>%
  distinct() %>%
  dplyr::filter()


all_log2fc_gpc3_ls <-
  list(Others = de_df_gpc3) %>%
  map(rownames_to_column, "gene") %>%
  map(function(all_marker){
    log2fc <- all_marker$avg_log2FC
    names(log2fc) <- all_gene_df$ENTREZID[match(all_marker$gene, all_gene_df$SYMBOL)]
    log2fc <- sort(log2fc, decreasing = TRUE)
    log2fc[!is.na(names(log2fc))]
  })  


all_gsea_gpc3_go_full_ls <- 
  all_log2fc_gpc3_ls %>%
  map(function(x){
    clusterProfiler::gseGO(
      geneList = x,
      ont = "ALL",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
  })


all_gsea_gpc3_kegg_full_ls <-   
  all_log2fc_gpc3_ls %>%
  map(function(x){
    clusterProfiler::gseKEGG(
      geneList = x,
      organism = "hsa",
      keyType = "kegg",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1
    )
  })


all_gsea_gpc3_kegg_full_ls$Others@result %>% view


# export
all_gsea_gpc3_kegg_full_ls %>%
  map('result') %>%
  purrr::discard(function(x) length(x) < 1) %>%
  imap(function(x, y){x$cell <- y; x}) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v5, "GSEA - enrich table of KEGG.xlsx"),
    asTable = TRUE
  )

# export
all_gsea_gpc3_go_full_ls %>%
  map('result') %>%
  purrr::discard(function(x) length(x) < 1) %>%
  imap(function(x, y){x$cell <- y; x}) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v5, "GSEA - enrich table of GO.xlsx"),
    asTable = TRUE
  )

## GSEA: KEGG and GO - GPC3 - PR/PD_SD -----------
# 20251104 add

all_log2fc_gpc3_ls_pos <-
  list(GPC3 = all_marker_gpc3_pos) %>%
  map(rownames_to_column, "gene") %>%
  map(function(all_marker){
    log2fc <- all_marker$avg_log2FC
    names(log2fc) <- all_gene_df$ENTREZID[match(all_marker$gene, all_gene_df$SYMBOL)]
    log2fc <- sort(log2fc, decreasing = TRUE)
    log2fc[!is.na(names(log2fc))]
  })  


all_gsea_gpc3_go_full_ls_pos <- 
  all_log2fc_gpc3_ls_pos %>%
  map(function(x){
    clusterProfiler::gseGO(
      geneList = x,
      ont = "ALL",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
  })


all_gsea_gpc3_kegg_full_ls_pos <-   
  all_log2fc_gpc3_ls_pos %>%
  map(function(x){
    clusterProfiler::gseKEGG(
      geneList = x,
      organism = "hsa",
      keyType = "kegg",
      exponent = 1,
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 1
    )
  })



# export
all_gsea_gpc3_kegg_full_ls_pos %>%
  map('result') %>%
  purrr::discard(function(x) length(x) < 1) %>%
  imap(function(x, y){x$cell <- y; x}) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v5, "GSEA - enrich table of KEGG - GPC3.xlsx"),
    asTable = TRUE
  )

# export
all_gsea_gpc3_go_full_ls_pos %>%
  map('result') %>%
  purrr::discard(function(x) length(x) < 1) %>%
  imap(function(x, y){x$cell <- y; x}) %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v5, "GSEA - enrich table of GO - GPC3.xlsx"),
    asTable = TRUE
  )

## nn cell props --------------
res_nn_name_df_gpc3  <-
  do.call(rbind, res_nn_name_ls) %>%
  as.data.frame() %>%
  rownames_to_column('obj_id') %>%
  left_join(
    all_data_meta_for_cn_analysis [c('cell_id_sample', 'celltype', 'project', 'group', 'annotation', 'annotation_full', 'obj_id')],
    by = "obj_id"
  ) %>%
    left_join(
      all_data_epi2[[c('cell_id_sample', 'AUC_GPC3', 'AUC_GPC3_Group')]],
      by = "cell_id_sample"
    )

res_nn_name_df_gpc3_full <-
  res_nn_name_df_gpc3 %>%
  pivot_longer(
    V1:V100,
    names_to = "NN_index",
    values_to = "NN_cells"
  )


# plot
# summarise by cluster
p_nn_props_col_gpc3 <-
  res_nn_name_df_gpc3_full %>%
  group_by(AUC_GPC3_Group, NN_cells) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = AUC_GPC3_Group, y = freq)) +
  geom_col(aes(fill = NN_cells)) +
  
  scale_fill_manual(values = col_plan_2) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )



p_nn_props_col_gpc3_group <-
  res_nn_name_df_gpc3_full %>%
  group_by(group, NN_cells) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = group, y = freq)) +
  geom_col(aes(fill = NN_cells)) +
  
  scale_fill_manual(values = col_plan_2) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )


p_nn_props_col_gpc3_group_double <-
  res_nn_name_df_gpc3_full %>%
  group_by(group, AUC_GPC3_Group, NN_cells) %>%
  summarise(count = n()) %>%
  group_by(group, AUC_GPC3_Group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = AUC_GPC3_Group, y = freq)) +
  geom_col(aes(fill = NN_cells)) +
  facet_wrap(~ group) +
  
  scale_fill_manual(values = col_plan_2) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )




ggsave(
  p_nn_props_col_gpc3,
  filename = file.path(out_dir_v5, "nn cells - bar plot of gpc3.pdf"), 
  height = 6,
  width  = 6
)

ggsave(
  p_nn_props_col_gpc3_group,
  filename = file.path(out_dir_v5, "nn cells - bar plot of group.pdf"), 
  height = 6,
  width  = 6
)

ggsave(
  p_nn_props_col_gpc3_group_double,
  filename = file.path(out_dir_v5, "nn cells - bar plot of gpc3 and group.pdf"), 
  height = 6,
  width  = 10
)



## 20251113 combine all immune cells to immune
p_nn_props_col_gpc3_group_double_combined_immune <-
  res_nn_name_df_gpc3_full %>%
  
  mutate(NN_cells = ifelse(
    NN_cells %in% c("Neutrophil", "Macrophage", "T", "B/Plasma cell", "Kupffer cell"),
    "Immune cell",
    NN_cells
  )) %>% 
  
  group_by(group, AUC_GPC3_Group, NN_cells) %>%
  summarise(count = n()) %>%
  group_by(group, AUC_GPC3_Group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = AUC_GPC3_Group, y = freq)) +
  geom_col(aes(fill = NN_cells)) +
  facet_wrap(~ group) +
  
  scale_fill_manual(values = col_plan_2) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )

ggsave(
  p_nn_props_col_gpc3_group_double_combined_immune,
  filename = file.path(out_dir_v5, "nn cells - bar plot of gpc3 and group - combined immune cells.pdf"), 
  height = 6,
  width  = 10
)


# DotPlot(
#   all_data_epi2,
#   group.by = "AUC_GPC3_Group",
#   features = gpc_de_gene_ls$GPC3_GROUP$gene
# ) +
#   coord_flip()
# 
# all_data_epi2 %>%
#   FetchData(
#     vars = gpc_de_gene_ls$GPC3_GROUP$gene %>% c("AUC_GPC3_Group")
#   ) %>%
#   group_by(AUC_GPC3_Group) %>%
#   summarise(across(everything(), mean)) %>%
#   column_to_rownames("AUC_GPC3_Group") %>%
#   apply(2, function(x) x[1] > x[2]) %>%
#   table


# all_data_cn %>%
#   DimPlot(
#     group.by = "RNA_snn_res.0.1",
#     reduction = "umap"
#   )
# 
# all_data_cn[[]] %>%
#   group_by(project, group, RNA_snn_res.0.1) %>% 
#   summarise(count = n()) %>%
#   group_by(RNA_snn_res.0.1) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   mutate(project = str_c(project, group, sep = "_")) %>% 
#   
#   ggplot(aes(x = RNA_snn_res.0.1, y = freq, fill = project)) +
#   geom_col() +
#   ggsci::scale_fill_igv()
# 
#   
# all_data_cn[[]] %>%
#   group_by(project, group, RNA_snn_res.0.1) %>% 
#   summarise(count = n()) %>%
#   group_by(project, group) %>%
#   mutate(freq = count/sum(count)) %>%
#   
#   mutate(project = str_c(project, group, sep = "_")) %>% 
#   
#   ggplot(aes(x =  project, y = freq, fill = RNA_snn_res.0.1)) +
#   geom_col() +
#   ggsci::scale_fill_igv()
  
## progeny -----------

### heatmap -----
all_data_epi2_progeny <-
  progeny::progeny(
    all_data_epi2,
    scale = FALSE,
    assay_name = "RNA",
    top = 100,
    perm = 1,
    return_assay = TRUE
  )

all_data_epi2_progeny <- all_data_epi2_progeny %>% ScaleData(assay = 'progeny')

df_progeny_gpc3 <-
  all_data_epi2_progeny %>%
  FetchData(
    vars = c(
      rownames(all_data_epi2_progeny[['progeny']]),
      'AUC_GPC3_Group', 
      'project',
      'group'
    )
  )


df_progeny_gpc3_tidy <- 
  df_progeny_gpc3 %>%
  group_by(project, group, AUC_GPC3_Group) %>%
  summarise(across(everything(), mean)) %>%
  dplyr::filter(group != "pre") 


p_heat_progeny_gpc3_ls  <-
  df_progeny_gpc3_tidy %>%
  group_by(AUC_GPC3_Group) %>%
  nest() %>%
  deframe() %>%
  # .[1] %>%
  imap(function(dat, name){ # browser()
    dat2 <- dat %>% arrange(group) %>% column_to_rownames("project") 
    colnames(dat2) <- dat2 %>% colnames() %>% str_remove("^progeny_")
    
    expr <- dat2[, ! colnames(dat2) %in% c("group")]  %>% scale  
    ht <- 
      ComplexHeatmap::Heatmap(
        expr,
        name = " ",
        col = colorRampPalette(rev(RColorBrewer::brewer.pal(Inf, "RdBu")))(50),
        left_annotation = ComplexHeatmap::rowAnnotation(
          df = dat2['group'],
          col = list(group = c(pre = "white", PR = "blue", "PD_SD" = "red"))
        ),
        cluster_columns = FALSE,
        row_split = dat2$group,
        column_title = name
      )
    ht
  })

# export
pdf(
  file = file.path(out_dir_v5, 'progeny pathway heatmap.pdf'),
  height = 4,
  width  = 8
)
p_heat_progeny_gpc3_ls %>% walk(print)
dev.off()  


### boxplot ---------
p_box_progeny_gpc3 <-
  df_progeny_gpc3_tidy %>%
  # dplyr::filter(AUC_GPC3_Group == "Others") %>%
  
  pivot_longer(
    starts_with("progeny"),
    names_to = "pathway",
    values_to = "val"
  ) %>%
  
  ggplot(aes(x = pathway, y = val, fill = group)) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75), shape = 21) +
  RotatedAxis() +
  facet_wrap(~ AUC_GPC3_Group) + 
  scale_y_log10() +
  ggpubr::geom_pwc(tip.length = 0)

# export
pdf(
  file = file.path(out_dir_v5, 'progeny pathway boxplot.pdf'),
  height = 4,
  width  = 12
)
p_box_progeny_gpc3 %>% print
dev.off()  



# 20251104 需求----------------
res_nn_ls <- readRDS("20251022 epi nn analysis/res_nn_ls.rds")
all_data_epi <- readRDS("all_data_epi.rds")


## tumor C0中GPC3/others占比-----------
p_nn_cluster_gpc3 <-
  all_data_epi[[]] %>%
  group_by(annotation_epi, AUC_GPC3_Group) %>%
  summarise(count = n()) %>%
  group_by(annotation_epi) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation_epi, y = freq, fill = AUC_GPC3_Group)) +
  geom_col() +
  ggsci::scale_fill_aaas()

ggsave(
  p_nn_cluster_gpc3,
  filename = file.path(out_dir_v5, str_glue("bar plot of gpc3 and others cells in each tumor nn cluster.pdf")),
  height = 6,
  width  = 8
)

## nn cells of tumor C0 ------
# 比较tumor C0周围的免疫细胞在pr pd中的差异情况

res_nn_name_ls <-
  res_nn_ls %>%
  map(function(x){ # browser()
    x$nn.idx %>% dim()
    
    x$nn.idx[1:20, 1:5]
    all_idx <- c(x$nn.idx) + 1
    all_celltype <- c("Others", as.character(x$dat$annotation_full))
    
    all_idx2celltype <- all_celltype[all_idx]
    all_idx2celltype_mat <-  
      matrix(
        all_idx2celltype,
        ncol = 100,
        byrow = FALSE
      )
    rownames(all_idx2celltype_mat) <- x$query_idx
    all_idx2celltype_mat
  })

res_nn_name_df <- 
  do.call(rbind, res_nn_name_ls) %>%
  as.data.frame() %>%
  rownames_to_column('obj_id') %>%
  left_join(
    all_data_meta_for_cn_analysis [c('celltype', 'project', 'group', 'annotation', 'annotation_full', 'obj_id')],
    by = "obj_id"
  ) %>%
  left_join(
    all_data_cn[[c('RNA_snn_res.0.1', 'obj_id')]],
    by = "obj_id"
  ) %>%
  mutate(tumor = str_c('Tumor C' , as.character(RNA_snn_res.0.1), sep = ""))

res_nn_name_df_full <-
  res_nn_name_df %>%
  pivot_longer(
    V1:V100,
    names_to = "NN_index",
    values_to = "NN_cells"
  ) 



p_nn_props_pie_tumor_grp_tumor <-
  res_nn_name_df_full %>%
  group_by(group, tumor, NN_cells) %>%
  summarise(count = n()) %>%
  group_by(group, tumor) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = 1, y = freq, fill = NN_cells)) +
  geom_col() +
  
  facet_grid(group ~ tumor) +
  coord_polar(theta = "y") +
  theme_classic() +
  scale_fill_manual(values = col_plan_2) +
  scale_y_continuous(labels = scales::percent) + 
  theme(
    strip.background = element_blank(),
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(
    x = NULL,
    y = NULL
  )

ggsave(
  p_nn_props_pie_tumor_grp_tumor,
  filename = file.path(out_dir_v5, str_glue("pie plot of nn cells - group and tumor cluster.pdf")),
  height = 8,
  width  = 10
)

## nn cells of PR / PD -----------------
all_grps_epi_grp_tumor <- 
  all_data_epi[[]] %>%
  group_by(group, annotation_epi) %>%
  nest() %>%
  mutate(
    grp = str_c(group, annotation_epi, sep = "---")
  ) %>%
  ungroup() %>%
  dplyr::select(grp, data) %>%
  deframe()


# tumor_c0_pr <- 
#   all_data_epi[[]]  %>%
#   dplyr::filter(
#     annotation_epi == "Tumor C0",
#     group == "PR"
#   ) %>%
#   pull('cell_id_sample')


res_nn_cell_id_ls <- res_nn_ls %>% # .[3] %>%
  imap(function(x, sample){ # browser()
    dat <- x$dat 
    
    query_cell_id <- unlist(dat[match(x$query_idx, dat$cell_id), 'cell_id_sample'])
    dat_cell_id   <- dat$cell_id_sample
    
    
    tumor_c0_pr <- all_grps_epi_grp_tumor$`PR---Tumor C0`$cell_id_sample
    tumor_c0_pd <- all_grps_epi_grp_tumor$`PD_SD---Tumor C0`$cell_id_sample
    
    table( query_cell_id %in%  tumor_c0_pr)
    table( query_cell_id %in%  tumor_c0_pd)
    
    tumor_c0_pr_cel_id <- tumor_c0_pr[tumor_c0_pr %in% query_cell_id]
    tumor_c0_pd_cel_id <- tumor_c0_pd[tumor_c0_pd %in% query_cell_id]
    
    nn_mat <- matrix(
      dat_cell_id[x$nn.dist],
      nrow = length(x$query_idx),
      ncol = 100,
      byrow = FALSE
    )
    rownames(nn_mat) <- query_cell_id
    
    list(
      nn_mat = nn_mat,
      pr_nn_mat = nn_mat[tumor_c0_pr_cel_id, ],
      pd_nn_mat = nn_mat[tumor_c0_pd_cel_id, ]
    )
  })


pr_nn_mat_cell_id <- 
  res_nn_cell_id_ls %>%
  map('pr_nn_mat') %>%
  # map_lgl(function(x) length(x) < 1) %>%
  do.call(rbind, .) %>%
  c() %>%
  unique()

pd_nn_mat_cell_id <- 
  res_nn_cell_id_ls %>%
  map('pd_nn_mat') %>%
  # map_lgl(function(x) length(x) < 1) %>%
  do.call(rbind, .) %>%
  c() %>%
  unique()

all_data_for_nn <- 
  all_data[, all_data$cell_id_sample %in% c(
    pr_nn_mat_cell_id,
    pd_nn_mat_cell_id
  )]

# all_data_for_nn$group_nn <- 
#   ifelse(
#     all_data_for_nn$cell_id_sample %in% pr_nn_mat_cell_id,
#     "PR",
#     "PD_SD"
#   )
# all( all_data_for_nn$group == all_data_for_nn$group_nn )

Idents(all_data_for_nn) <- "group"
all_marker_nn <- 
  FindMarkers(
    all_data_for_nn,
    ident.1 = "PR",
    ident.2 = "PD_SD",
    logfc.threshold = 0,
    min.pct = 0
  )

# all_marker_nn
write.csv(
  all_marker_nn,
  file = file.path(out_dir_v5, "de analysis of nn cells in Tumor C0 between PR and PD_SD.csv"),
  row.names = TRUE
)

all_de_genes_for_enrichmenet <-
  all_marker_nn %>% 
  rownames_to_column('gene') %>%
  mutate(cluster = 'nn') %>%
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map(function(x){
    up_gene <- 
      x %>%
      dplyr::filter(
        p_val_adj < 0.05,
        # pct.1 - pct.2 > 0.1,
        avg_log2FC > 0.25
      ) %>%
      .$gene
    
    down_gene <- 
      x %>% 
      dplyr::filter(
        p_val_adj < 0.05,
        # pct.1 - pct.2 > 0.1,
        avg_log2FC < -0.25
      ) %>%
      .$gene
    
    # return
    list(
      PR = up_gene,
      PD_SD = down_gene
    )
  })  %>%
  flatDoubleList()

all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)

### enrichemt analysis -----
gene_db <- 
  clusterProfiler::bitr(
    # unique(unlist(all_de_genes_id_for_enrichmenet)),
    rownames(all_data_for_nn),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = "org.Hs.eg.db"
  ) 


dir_enrich <- file.path(out_dir_v5, "nn enrichment analysis")
dir.create(dir_enrich)

all_de_genes_id_for_enrichmenet <-
  all_de_genes_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::bitr(
      x,
      OrgDb = "org.Hs.eg.db",
      fromType = "SYMBOL",
      toType = "ENTREZID"
    ) %>%
      .$ENTREZID
  }) %>%
  discard(function(x) length(x) < 1)


#### KEGG ------------
enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      x,
      organism = "hsa", 
      keyType = "ncbi-geneid", 
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      universe = gene_db$ENTREZID,
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2,
      use_internal_data =  FALSE
    )
  })

enrich_ls <-
  enrich_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    # gene symbol
    xx2 <-
      xx %>%
      map(function(g){ # browser()
        gene_name_vec <- gene_db$SYMBOL[match(g, gene_db$ENTREZID)] 
        gene_name <- str_c(str_replace_na(gene_name_vec), collapse = ",")
        
        # gene_loc_vec <- gene_description$location[match(gene_name_vec, gene_description$gene_name)]
        # gene_loc  <- str_c(str_replace_na(gene_loc_vec), collapse = ",")
        
        list(
          gene_name
          # gene_loc
        )
      })
    
    x@result$geneSymbol <- xx2 %>% map(1) %>% unlist(use.names = FALSE)
    # x@result$geneLocation <- xx2 %>% map(2) %>% unlist(use.names = FALSE)
    x
  })

# saveRDS(enrich_ls, file.path(dir_enrich, 'enrich_ls.rds'))

# expot enrichment table
enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(dir_enrich, str_glue("KEGG - enriment table - {name}.csv")),
      row.names = FALSE
    )
  })

# dotplot
enrich_dot_ls <-
  enrich_ls %>%
  imap(function(x, name){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export dotplot
enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(dir_enrich, str_glue("dotplot - {name}.pdf")),
      height = 6,
      width = 8
    )
  })

#### GOBP -------
gobp_enrich_ls <-
  all_de_genes_id_for_enrichmenet %>%
  map(function(x){
    clusterProfiler::enrichGO(
      x,
      OrgDb = 'org.Hs.eg.db',
      ont = "BP",
      pvalueCutoff = 0.05, 
      pAdjustMethod = "BH", 
      minGSSize = 10,
      maxGSSize = 500, 
      qvalueCutoff = 0.2
    )
  })

gobp_enrich_ls <-
  gobp_enrich_ls %>%
  map(function(x){ # browser()
    xx <- x@result$geneID %>% str_split("/")
    
    # gene symbol
    xx2 <-
      xx %>%
      map(function(g){ # browser()
        gene_name_vec <- gene_db$SYMBOL[match(g, gene_db$ENTREZID)] 
        gene_name <- str_c(str_replace_na(gene_name_vec), collapse = ",")

        # gene_loc_vec <- gene_description$location[match(gene_name_vec, gene_description$gene_name)]
        # gene_loc  <- str_c(str_replace_na(gene_loc_vec), collapse = ",")
        
        list(
          gene_name
          # gene_loc
        )
      })
    
    x@result$geneSymbol <- xx2 %>% map(1) %>% unlist(use.names = FALSE)
    # x@result$geneLocation <- xx2 %>% map(2) %>% unlist(use.names = FALSE)
    x
  })
# saveRDS(gobp_enrich_ls, file.path(dir_enrich, 'gobp_enrich_ls.rds'))

# expot enrichment table
gobp_enrich_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    
    write.csv(
      x@result,
      file = file.path(dir_enrich, str_glue("GOBP - enriment table - {name}.csv")),
      row.names = FALSE
    )
  })

# dotplot
gobp_enrich_dot_ls <-
  gobp_enrich_ls %>%
  imap(function(x, name){
    if(nrow(as.data.frame(x)) < 1){
      return(NULL)
    }
    enrichplot::dotplot(
      x,
      showCategory = 10,
      label_format = 50
    ) +
      ggtitle(name)
  })

# export dotplot
gobp_enrich_dot_ls %>%
  discard(function(x) length(x) < 1) %>%
  iwalk(function(x, name){
    ggsave(
      x,
      filename = file.path(dir_enrich, str_glue("BOBP - dotplot - {name}.pdf")),
      height = 6,
      width = 8
    )
  })










