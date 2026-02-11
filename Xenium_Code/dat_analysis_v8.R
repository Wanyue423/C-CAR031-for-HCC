library(tidyverse)
library(Seurat)

# 需求： 20251120需求.docx

all_data2_cd45 <- readRDS('all_data2_cd45.rds')

all_data2 <- readRDS('all_data2.rds')
all_data2$project <- all_data2$sample

# all_data2_cd45 <- readRDS('all_data2_cd45.rds')
# all_data2_t <- all_data2_cd45[, all_data2_cd45$annotation == "T"]


out_dir_v8 <- "20251120 dat analysis"
dir.create(out_dir_v8)


# 1. dotplot of selected pathway --------------
select_path_ls <- list(
  Immunosuppression = c('ICOS', 'LGALS9', 'CD276', 'CD160', 'ARG1', 'LAG3', 'VTCN1', 'CD27', 'CD28'),
  Chemoattractant   = c('IFNG', 'IL7', 'IL2RA', 'TNFSF8', 'CCL20'),
  APC               = c('CD40', 'CD80', 'CD86'),
  'ECM remodeling'  = c('SERPINE1', 'THBS1', 'PDGFA')
)
select_path_ls <- 
  select_path_ls %>%
  enframe() %>%
  unnest(value) %>%
  deframe()

# macrophage
sc_mac <- all_data2[, all_data2$annotation == "Macrophage" & all_data2$group != "pre"]
p_dot_mac_select <-
  DotPlot(
    sc_mac,
    group.by = 'group',
    features = select_path_ls
  ) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(Inf, "RdBu") %>% {.[-c(1, length(.))]} %>% rev() ) +
  scale_x_discrete(expand = c(0, 1)) +
  theme(
    strip.clip = 'off',
    panel.spacing = unit(0, 'mm'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL
  )

ggsave(
  p_dot_mac_select,
  filename = file.path(out_dir_v8, '1 dotplot of selected gene - macrophage.pdf'),
  height = 4,
  width = 12
)

# neutrophil
sc_neu <- all_data2[, all_data2$annotation == "Neutrophil" & all_data2$group != "pre"]
p_dot_neu_select <-
  DotPlot(
    sc_neu,
    group.by = 'group',
    features = select_path_ls
  ) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(Inf, "RdBu") %>% {.[-c(1, length(.))]} %>%  rev()) +
  scale_x_discrete(expand = c(0, 1)) +
  theme(
    strip.clip = 'off',
    panel.spacing = unit(0, 'mm'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL
  )

ggsave(
  p_dot_neu_select,
  filename = file.path(out_dir_v8, '1 dotplot of selected gene - neutrophil.pdf'),
  height = 4,
  width = 12
)



# 2. macrophage gsea result --------- 
all_markers_major_cd45_ls <-
  all_data2_cd45$annotation %>%
  unique() %>%
  setNames(., .) %>%
  map(function(x){
    Idents(all_data2_cd45) <- "annotation"
    all_markers <- 
      FindMarkers(
        all_data2_cd45,
        ident.1 = "PR",
        ident.2 = "PD_SD",
        group.by = "group",
        subset.ident = x,
        logfc.threshold = 0
      )
    all_markers$gene <- rownames(all_markers)
    all_markers$cluster <- x
    all_markers
  }) 


all_gene_df <- 
  clusterProfiler::bitr(
    rownames(all_data2_cd45),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ) %>%
  distinct() %>%
  dplyr::filter()


all_log2fc_cd45_major_ls <-
  all_markers_major_cd45_ls %>%
  map(function(all_marker){
    log2fc <- all_marker$avg_log2FC
    names(log2fc) <- all_gene_df$ENTREZID[match(all_marker$gene, all_gene_df$SYMBOL)]
    log2fc <- sort(log2fc, decreasing = TRUE)
    log2fc[!is.na(names(log2fc))]
  })



all_gsea_cd45_go_mac_ls <- 
  all_log2fc_cd45_major_ls['Macrophage'] %>%
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

# dotplot
gsea_go_mac <- all_gsea_cd45_go_mac_ls$Macrophage
gsea_go_mac@result <- 
  gsea_go_mac@result  %>%
  dplyr::filter(
    ONTOLOGY == 'BP',
    p.adjust < 0.05
  ) 

p_dot_mac <-
  gsea_go_mac %>%
  enrichplot::dotplot(
    showCategory = 15,
    label_format = 50,
    # orderBy = 'p.adjust'
  ) +
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(Inf, "RdBu") %>% {.[-c(1, length(.))]} %>%  rev()) +
  theme(
    panel.grid = element_blank()
  )

# export
ggsave(
  p_dot_mac,
  filename = file.path(out_dir_v8, '2 dotplot of signicant GSEA GO BP pathway - macrphage.pdf'),
  height = 6,
  width = 8
)


# barplot
p_gsea_barplot_mac <-
  gsea_go_mac@result %>%
  mutate(group = ifelse(NES < 0, "PD enrichment", "PR enrichment")) %>% 
  arrange(group, abs(NES)) %>%
  mutate(Description = fct_inorder(Description)) %>%
  # mutate(
  #   pvalue = -log10(pvalue),
  #   pvalue = pmax(pvalue, 0.1),
  #   pvalue = ifelse(group == 'PD', -pvalue, pvalue)
  # ) %>%
  ggplot(aes(x = NES, y = Description, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(
    data = ~ filter(., group == 'PR enrichment'),
    aes(x = 0, label = Description),
    hjust = 0,
    size = 3
  ) +
  geom_text(
    data = ~ filter(., group == 'PD enrichment'),
    aes(x = 0, label = Description),
    hjust = 1,
    size = 3
  ) +
  geom_vline(xintercept = 0) +
  ggsci::scale_fill_bmj() +
  scale_x_continuous(
    limits = function(x) c(median(x) + 1 * c(-1, 1.5) * diff(range(x)) ),
    labels = function(x) abs(x)
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.9, 0.2),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 15)
  ) +
  coord_cartesian(clip = 'off') +
  labs(
    x = 'Normalized Enrichment Score'
  )

# export
ggsave(
  p_gsea_barplot_mac,
  filename = file.path(out_dir_v8, '2 barplot of signicant GSEA GO BP pathway - macrphage.pdf'),
  height = 6,
  width = 8
)

# 20251211 change color 
p_gsea_barplot_mac2 <- 
  p_gsea_barplot_mac +
  scale_fill_manual(values = c(
    "PD enrichment" = "#2166AC",
    "PR enrichment" = "#B2182B"
  ))

ggsave(
  p_gsea_barplot_mac2,
  filename = file.path(out_dir_v8, '2 barplot of signicant GSEA GO BP pathway - macrphage - v2.pdf'),
  height = 6,
  width = 8
)

# 20251211 change color - v2
p_gsea_barplot_mac3 <-
  gsea_go_mac@result %>% 
  ggplot(aes(x = NES, y = reorder(Description, NES), fill = NES)) + 
  geom_col() + 
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu"))) + 
  theme( 
    legend.position = "none" 
  ) + 
  labs( 
    y = NULL 
  ) 

ggsave(
  p_gsea_barplot_mac3,
  filename = file.path(out_dir_v8, '2 barplot of signicant GSEA GO BP pathway - macrphage - v3.pdf'),
  height = 6,
  width = 8
)

# 3. nn cluster of epithelial ------------
# all_data_cn <- readRDS("20251022 epi nn analysis/all_data_cn.rds")
all_data_cn2 <- all_data_cn[, all_data_cn$RNA_snn_res.0.1 != '5']

pdf(
  file = file.path(out_dir_v8, "3 cn cluster - umap plot of cluster.pdf"),
  height = 6,
  width  = 8
)
DimPlot(all_data_cn2, group.by = "RNA_snn_res.0.1", label = TRUE) + NoLegend() 
dev.off()


pdf(
  file = file.path(out_dir_v8, "3 cn cluster - celltype(marker) heatmap of cluster.pdf"),
  height = 5,
  width  = 10
)
AverageExpression(
  all_data_cn2,
  group.by = "RNA_snn_res.0.1"
) %>%
  .$RNA %>%
  as.matrix() %>% t %>%
  pheatmap::pheatmap(
    scale = "column",
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    border_color = NA,
    labels_row = str_c("Tumor C", seq.int(nrow(.)) - 1, sep = "")
  )
dev.off()

pdf(
  file = file.path(out_dir_v8, "3 cn cluster - tumor cluster prop of group - add pre.pdf"),
  height = 6,
  width  = 4
)
all_data_cn2[[]] %>% 
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


# 4. barplot of GSEA tumor C0 -------------
gsea_tumorc0_pd <- read.csv("20251022 epi nn analysis/nn enrichment analysis/KEGG - enriment table - nn---PD_SD.csv")
gsea_tumorc0_pr <- read.csv("20251022 epi nn analysis/nn enrichment analysis/KEGG - enriment table - nn---PR.csv")


gsea_tumorc0_pd_filter <-
  gsea_tumorc0_pd %>%
  dplyr::filter(ID %in% c(
    'hsa01100',
    'hsa04610',
    'hsa00982',
    'hsa04010',
    'hsa04151',
    'hsa04350'
  ))

gsea_tumorc0_pr_filter <-
  gsea_tumorc0_pr %>%
  dplyr::filter(ID %in% c(
    'hsa04670',
    'hsa04064',
    'hsa04668',
    'hsa04060'
  ))


gsea_tumorc0 <-
  list(
    PD = gsea_tumorc0_pd_filter,
    PR = gsea_tumorc0_pr_filter
  ) %>%
  enframe('group') %>%
  unnest(value)


p_gsea_tumorc0 <-
  gsea_tumorc0 %>%
  arrange(group, desc(pvalue)) %>%
  mutate(Description = fct_inorder(Description)) %>%
  mutate(
    pvalue = -log10(pvalue),
    pvalue = pmax(pvalue, 0.1),
    pvalue = ifelse(group == 'PD', -pvalue, pvalue)
  ) %>%
  ggplot(aes(x = pvalue, y = Description, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(
    data = ~ filter(., group == 'PR'),
    aes(x = 0, label = Description),
    hjust = 0,
    size = 3
  ) +
  geom_text(
    data = ~ filter(., group == 'PD'),
    aes(x = 0, label = Description),
    hjust = 1,
    size = 3
  ) +
  geom_vline(xintercept = 0) +
  ggsci::scale_fill_bmj() +
  scale_x_continuous(
    limits = function(x) c(median(x) + 1 * c(-1, 1.5) * diff(range(x)) ),
    labels = function(x) abs(x)
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.9, 0.2),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 15)
  ) +
  coord_cartesian(clip = 'off') +
  labs(
    x = expression( log[10] * '(' * 'P' * ')' )
  )


# export
ggsave(
  p_gsea_tumorc0,
  filename = file.path(out_dir_v8, '4 barplot of GSEA tumor C0.pdf'),
  height = 6,
  width = 6
)



# 5. barplot of GSEA GPC3 ------
gsea_gpc3 <- readxl::read_excel('20251022 epi nn analysis/GSEA - enrich table of KEGG - GPC3.xlsx')

gsea_gpc3_filter <-
  gsea_gpc3 %>%
  dplyr::filter(ID %in% c(
    'hsa04064',
    'hsa04370',
    'hsa04668',
    'hsa04215',
    
    'hsa04010',
    'hsa04150',
    'hsa04310',
    'hsa04151'
  ))
p_gsea_gpc3 <-
  gsea_gpc3_filter %>%
  mutate(group = ifelse(NES < 0, "Neg", "Pos")) %>% 
  arrange(group, abs(NES)) %>%
  mutate(Description = fct_inorder(Description)) %>%
  # mutate(
  #   pvalue = -log10(pvalue),
  #   pvalue = pmax(pvalue, 0.1),
  #   pvalue = ifelse(group == 'PD', -pvalue, pvalue)
  # ) %>%
  ggplot(aes(x = NES, y = Description, fill = group)) +
  geom_col(width = 0.75) +
  geom_text(
    data = ~ filter(., group == 'Pos'),
    aes(x = 0, label = Description),
    hjust = 0,
    size = 3
  ) +
  geom_text(
    data = ~ filter(., group == 'Neg'),
    aes(x = 0, label = Description),
    hjust = 1,
    size = 3
  ) +
  geom_vline(xintercept = 0) +
  ggsci::scale_fill_bmj() +
  scale_x_continuous(
    limits = function(x) c(median(x) + 1 * c(-1, 1.5) * diff(range(x)) ),
    labels = function(x) abs(x)
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = c(0.9, 0.2),
    axis.text.x = element_text(size = 12),
    legend.text = element_text(size = 10),
    axis.title.x = element_text(size = 15)
  ) +
  coord_cartesian(clip = 'off') +
  labs(
    x = 'Normalized Enrichment Score'
  )


# export
ggsave(
  p_gsea_gpc3,
  filename = file.path(out_dir_v8, '5 barplot of GSEA tumor GPC3.pdf'),
  height = 6,
  width = 6
)





