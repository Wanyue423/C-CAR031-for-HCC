library(tidyverse)
library(Seurat)

# all_data <- readRDS("report_all_samples/output/data/adata_celltype.rds")

all_data2 <- readRDS('all_data2.rds')
all_data2$project <- all_data2$sample

out_dir_v10 <- "20251209 dat analysis"
dir.create(out_dir_v10)



# 1 marker dotplot ------------

# all_data$res_1.0 %>% unique


# DimPlot(
#   all_data,
#   group.by = 'res_1.0'
# )
# 
# 
# DotPlot(
#   all_data,
#   group.by = "res_1.0",
#   features = c('UBD')
# )




marker_ls <- c(
  "B cell" = "CD19", 
  "B cell" = "CD79A",
  "B cell" = "MS4A1",

  "Plasma" = "TNFRSF17", 
  "Plasma" = "MZB1", 
  
  "Epithelial" = "EPCAM",
  # "Epithelial" = "KRT20", 
  "Cholangiocyte" = "KRT7", 
  "Cholangiocyte" = "KRT19", 
  "Cholangiocyte" = "SOX9",
  
  "Endothelial" = "PECAM1", 
  
  
  
  
  "Hepatic stellate cells" = 'ACTA2',
  "Hepatic stellate cells" = "PDGFRA", 
  "Hepatic stellate cells" = "CCDC80", 
  
  # "Smooth muscle cells" = "PDGFRB", 
  # "Smooth muscle cells" = "PLN", 
  # "Smooth muscle cells" = "MYLK", 
  # "Smooth muscle cells" = "TAGLN",
  
  
  
  
  "Hepatocyte" = "ALB", 
  "Hepatocyte" = "CYP3A4", 
  "Hepatocyte" = "CYP2E1", 
  "Hepatocyte" = "APOB", 
  "Hepatocyte" = "APOE",
  "Hepatocyte" = "HP",
  
  "Macrophage" = "CD68", 
  "Macrophage" = "CD163", 
  
  # "Kupffer cell" = "MAFB", 
  "Kupffer cell" = "VSIG4",
  
  
  "Liver sinusoidal\nendothelial cells" = "CLEC4G", 
  "Liver sinusoidal\nendothelial cells" = "CLEC4M",
  
  
  
  'Neutrophi' = "CSF3R", 
  'Neutrophi' = "FCGR3B", 
  
  "Pericyte" = "RGS5",
  
  
  "Proliferating" = "MKI67", 
  "Proliferating" = "TOP2A",
  
  
  "T cell" = "CD3E",
  "T cell" = "CD4",
  "T cell" = "FOXP3", 
  "T cell" = "CD8A", 
  "T cell" = "GZMB",
  
  
  # "NK" = "NCAM1", 
  
  "Tumor" = "UBD"
)



# DotPlot(
#   all_data,
#   group.by = "res_1.0",
#   features = marker_ls
# )




p_dot_marker <-
  DotPlot(
    all_data2,
    group.by = "annotation",
    features = marker_ls
  ) +
  RotatedAxis() +
  theme(
    strip.clip = "off",
    panel.spacing.x = unit(20, "mm")
  )


ggsave(
  p_dot_marker,
  filename = file.path(out_dir_v10, 'dotplot of marker.pdf'),
  height = 8,
  width = 28
)


# 2 ITH score of DEPTH analysis ---------
pre_grp <- c(
  PFS_long = "CJU-pre",
  PFS_long = "HYLO-pre",
  PFS_long = "XXLI-pre",
  PFS_long = "LBCH-pre",
  
  PFS_short = "XWRO-pre",
  PFS_short = "LZQI-pre",
  PFS_short = "LXZH-pre"
)

all_data_epi_pre <-
  all_data2[, 
    all_data2$sample %in% pre_grp &
    all_data2$annotation %in% c("Tumor cell", "Proliferating tumor cell")
  ]


# ITH
score_ith <-
  DEPTH::DEPTH(
    all_data_epi_pre@assays$RNA$data,
    match = data.frame(State = colnames(all_data_epi_pre), Identification = "Tumor")
  )

stopifnot(all(colnames(all_data_epi_pre) == score_ith[, 'Sample']))

all_data_epi_pre$ITH_score <- as.numeric(score_ith[, 'ITH score'])


# ITH score
p_ith_score <-
  all_data_epi_pre[[]] %>%
  dplyr::filter(sample %in% pre_grp) %>%
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>%
  ggplot(aes(x = PFS_group, y = ITH_score, fill = PFS_group)) +
  geom_violin(trim = F) +
  geom_boxplot(width = 0.15, fill = NA, outliers = F) +
  ggpubr::geom_pwc(tip.length = 0, bracket.nudge.y = 0.1) +
  ggpubr::theme_classic2()



ggsave(
  p_ith_score,
  filename = file.path(out_dir_v10, 'vlnplot of ITH score between PFS long and short.pdf'),
  height = 6,
  width = 8
)

    
# 3 enrichment analysis of GPC3 cells (PR vs PD) -----------

# file.path(out_dir_v5, "GPC3 cells - de analysis table of PR vs PD_SD in GPC3.csv")

de_marker_gpc3 <- 
  read.csv(
    file.path("20251022 epi nn analysis", "GPC3 cells - de analysis table of PR vs PD_SD in GPC3.csv"),
    row.names = 1
  )


de_gene_gpc3_ls <- 
  de_marker_gpc3 %>%
  rownames_to_column("gene") %>% 
  mutate(group = fct_recode(group, PR = 'Pos', PD = "Neg")) %>%
  dplyr::filter(group != "Non") %>%
  
  group_by(group) %>%
  nest() %>%
  deframe() %>%
  
  map('gene')




de_geneid_gpc3_ls <-
  de_gene_gpc3_ls %>%
  map(~{
    clusterProfiler::bitr(
      .x,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = "org.Hs.eg.db"
    ) %>%
      .$ENTREZID
  })


enrich_gpc3_ls <- 
  de_geneid_gpc3_ls %>%
  map(function(x){
    clusterProfiler::enrichGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

# enrichment 
enrich_gpc3_ls <-
  enrich_gpc3_ls %>%
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
p_enrich_gpc3_ls <-
  enrich_gpc3_ls %>%
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


# export
pdf(
  file.path(out_dir_v10, "enrichment analysis of gpc3 cell de gene of PR vs PD.pdf"), 
  height = 6.5, 
  width  = 8
)
p_enrich_gpc3_ls %>%  walk(print)
dev.off()

enrich_gpc3_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file.path(out_dir_v10, "enrichment analysis of gpc3 cell de gene of PR vs PD.xlsx"),
    asTable = TRUE
  )


## 3.2 KEGG enrichment ------------
# enrichment 
enrich_kegg_gpc3_ls <- 
  de_geneid_gpc3_ls %>%
  map(function(x){
    clusterProfiler::enrichKEGG(
      gene = x,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pvalueCutoff  = 1, # 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1 # 0.2
    )
  })

enrich_kegg_gpc3_ls <-
  enrich_kegg_gpc3_ls %>%
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
p_enrich_kegg_gpc3_ls <-
  enrich_kegg_gpc3_ls %>%
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


# export
pdf(
  file.path(out_dir_v10, "enrichment (KEGG) analysis of gpc3 cell de gene of PR vs PD.pdf"), 
  height = 6.5, 
  width  = 8
)
p_enrich_kegg_gpc3_ls %>%  walk(print)
dev.off()

enrich_kegg_gpc3_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file.path(out_dir_v10, "enrichment (KEGG) analysis of gpc3 cell de gene of PR vs PD.xlsx"),
    asTable = TRUE
  )


# 4. nn cluster of epithelial ------------

add_minicoord <- function(p, scale = 1, position = 0.1, size = 0.3, xlab = "UMAP1", ylab = "UMAP2"){
  require(cowplot)
  
  # 绘制坐标系
  arrow <- grid::arrow(angle = 20, type = "closed", length = unit(0.1, "npc"))
  
  umap_coord <- ggplot(tibble(group = c(xlab, ylab),
                              x = c(0, 0), xend = c(1, 0),
                              y = c(0, 0), yend = c(0, 1),
                              lx = c(0.4, -0.15), ly = c(-0.15, 0.4),
                              angle = c(0, 90))) +
    geom_segment(aes(x, y, xend = xend, yend = yend, group = group),
                 arrow = arrow, size = 1, lineend = "round", linejoin = "mitre") +
    geom_text(aes(lx, ly, label = group, angle = angle), size = 4) +
    theme_void() +
    coord_fixed(xlim = c(-0.3, 1), ylim = c(-0.3, 1))
  
  
  # 使用cowplot进行拼图
  if("patchwork" %in% class(p)){
    p_no_legend <- p & Seurat::NoAxes()
  }else{
    p_no_legend <- p + Seurat::NoAxes()
  }
  p_dim_full <-
    ggdraw() +
    draw_plot(
      p_no_legend,
      x = position, 
      y = position,
      width  = 1 - position,
      height = 1 - position,
      scale = scale
    ) +
    draw_plot(
      umap_coord,
      x = 0,
      y = 0,
      width = size,
      height = size
    )
  
  
  # 绘图即可
  p_dim_full
}

# all_data_cn <- readRDS("20251022 epi nn analysis/all_data_cn.rds")
# all_data_cn2 <- all_data_cn[, all_data_cn$RNA_snn_res.0.1 != '5']

p_umap_epi_nn <- 
  DimPlot(
    all_data_cn2, 
    group.by = "RNA_snn_res.0.1", 
    pt.size = 2,
    alpha = 0.1,
    label = TRUE
  ) +
NoLegend() 
 
set.seed(1234)
dat_mask <- 
  p_umap_epi_nn$data %>%
  slice_sample(n = 10000)

maskTable <- mascarade::generateMask(
  dims     = dat_mask[1:2], 
  clusters = dat_mask$RNA_snn_res.0.1,
  # expand = 1/200
)

maskTable2 <- 
  maskTable %>%
  group_by(cluster) %>%
  nest() %>%
  deframe() %>%
  map(function(dat){
    dat2 <- polyclip::polyoffset(
      list(x = dat$umap_1, y = dat$umap_2),
      delta = -0.1
    ) %>%
    .[[1]]
    names(dat2) <- c('umap_1', 'umap_2')
    as.data.frame(dat2)
  }) %>%
  enframe('cluster') %>%
  unnest(value)

p_umap_epi_nn2 <-
  p_umap_epi_nn +
  geom_path(
    data = maskTable2,
    aes(x = umap_1, y = umap_2, group = cluster, color = cluster),
    size = 1.5,
    linetype = 2
  ) + 
  ggsci::scale_color_bmj()



pdf(
  file = file.path(out_dir_v10, "cn cluster - umap plot of cluster - v3.pdf"),
  height = 6,
  width  = 8
)
add_minicoord(p_umap_epi_nn2)
dev.off()



# 5. umap plot of T subcluster ------------
# all_data_t <- readRDS('report_all_samples/output/data/adata_t.rds')
# all_data_t2 <- all_data_t[, str_detect(all_data_t$celltype, "like", negate = TRUE)]

/home/szt/miniforge3/envs/rpy/bin/python

import os
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rc_context

adata_t = sc.read('report_all_samples/output/data/adata_t.h5ad')

out_dir_v10 = "20251209 dat analysis"


# with rc_context({'figure.figsize':(6,6)}):
#   sc.pl.umap(
#     adata_t,
#     color='celltype', 
#     legend_loc='on data', 
#     legend_fontsize=10,
#     legend_fontoutline=2,
#     legend_fontweight='normal',
#     show=False
#   )
# plt.savefig(os.path.join(out_dir, 't_umap.pdf'))


with rc_context({'figure.figsize':(6,6)}):
  sc.pl.umap(
    adata_t,
    color='celltype',
    wspace = 0.1
  )
  plt.savefig(os.path.join(out_dir_v10, 't_umap.pdf'), dpi=250, bbox_inches='tight', facecolor='white')    





