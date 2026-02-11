library(tidyverse)
library(Seurat)

all_data <- readRDS("report_all_samples/output/data/adata_celltype.rds")

analysis_out_dir <- "out_dir_20250627"
dir.create(analysis_out_dir)

t_cells <- c(
  "Tregs",
  "Gamma delta T", 
  "Hepatocyte-like T cell", 
  "Macrophage-like T cell",
  "Fibroblast-like T cell",  
  "Cycling T cell",
  "CD4T",
  "CD8Tem"
)

# adata2
all_data2 <- all_data[, all_data$celltype != "low-quality cell"]
all_data2$annotation <- 
  ifelse(
    all_data2$celltype %in% t_cells,
    "T",
    as.character(all_data2$celltype)
  )

# saveRDS(all_data2, 'all_data2.rds')


cell_colors_used <- 
  all_data$celltype %>%
  levels() %>% 
  c(., "T") %>% 
  {structure(
    col_plan_1[seq_along(.)],
    names = .
  )}

grp_colors_used <- c(
  Responder = "#484848",
  'Non-responder' = "#EE8227"
)

#1. umap -------
p_dim_1 <- 
  DimPlot(
    all_data, 
    group.by = "celltype",
    label = TRUE,
    repel = TRUE,
  ) +
  NoLegend() +
  scale_color_manual(values = cell_colors_used)



p_dim_2 <- 
  DimPlot(
    all_data2, 
    group.by = "annotation",
    label = TRUE,
    repel = TRUE,
  ) +
  NoLegend() +
  scale_color_manual(values = cell_colors_used)

# export
ggsave(
  p_dim_1,
  filename = file.path(analysis_out_dir, "umap plot of annotation - all cells.pdf"),
  height = 7,
  width  = 8
)
ggsave(
  p_dim_2,
  filename = file.path(analysis_out_dir, "umap plot of annotation - remove low-quality cell and rename T cells.pdf"),
  height = 7,
  width  = 8
)


#2. vlnplot of gene expression ----------------

## tumor PR < PD_SD -----------



# p_vln_tumor_1 <-
#   all_data2[, all_data2$group != "pre" & all_data2$annotation == "Tumor cell"] %>%
#   VlnPlot(
#     features = c("TGFB1", "CXCL9"),
#     group.by = "group"
#   )
# 
# p_vln_tumor_2 <-
#   all_data2[, all_data2$group != "pre" & all_data2$annotation == "Proliferating tumor cell"] %>%
#   VlnPlot(
#     features = c("TGFB1", "CXCL9"),
#     group.by = "group"
#   )

dat_tgf_tumor <- 
  all_data2 %>%
  FetchData(vars = c("TGFB1", "CXCL9"), layer = "data") %>%
  cbind(all_data2[[c('celltype', 'group')]])


dat_tgf_tumor_tidy <-
  dat_tgf_tumor %>%
  dplyr::filter(group != "pre") %>% 
  dplyr::filter(celltype %in% c("Tumor cell", "Proliferating tumor cell")) %>%
  pivot_longer(
    all_of(c("TGFB1", "CXCL9")),
    names_to = "gene",
    values_to = "expr"
  )  %>%
  
  # 20250630 adjust group
  mutate(
    group = fct_recode(
      group,
      Responder = "PR",
      'Non-responder' = 'PD_SD'
    ) %>%
      fct_relevel(
        c('Responder', 'Non-responder')
      )
  )


dat_tgf_tumor_tidy2 <-
  dat_tgf_tumor_tidy %>%
  mutate(celltype = "Total tumor cells") %>% 
  rbind(dat_tgf_tumor_tidy) 


p_tgf_tumor_1 <- 
  dat_tgf_tumor_tidy2 %>%
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(celltype ~ gene) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  # scale_y_continuous(transform = "log10") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) +
  # 20250630 add
  scale_fill_manual(values = grp_colors_used)
 
p_tgf_tumor_2 <- 
  dat_tgf_tumor_tidy2 %>%
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(celltype ~ gene) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  scale_y_continuous(transform = "log1p") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) +
  # 20250630 add
  scale_fill_manual(values = grp_colors_used)

p_tgf_tumor_3 <-
  dat_tgf_tumor_tidy2 %>%
  
    dplyr::filter(expr != 0) %>%   
  
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(celltype ~ gene) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  # scale_y_continuous(transform = "log1p") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) +
  # 20250630 add
  scale_fill_manual(values = grp_colors_used)

p_tgf_tumor_4 <-
  dat_tgf_tumor_tidy2 %>%
  
  dplyr::filter(expr != 0) %>%   
  
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(celltype ~ gene) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  scale_y_continuous(transform = "log10") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) +
  # 20250630 add
  scale_fill_manual(values = grp_colors_used)

# export
ggsave(
  p_tgf_tumor_1,
  filename = file.path(analysis_out_dir, "vlnplot plot of TGFB1 and CXCL9 in tumor cells.pdf"),
  height = 10,
  width = 8
)
ggsave(
  p_tgf_tumor_2,
  filename = file.path(analysis_out_dir, "vlnplot plot of TGFB1 and CXCL9 in tumor cells - log transformed exression.pdf"),
  height = 10,
  width = 8
)
ggsave(
  p_tgf_tumor_3,
  filename = file.path(analysis_out_dir, "vlnplot plot of TGFB1 and CXCL9 in tumor cells - remove zero values.pdf"),
  height = 10,
  width = 8
)
ggsave(
  p_tgf_tumor_4,
  filename = file.path(analysis_out_dir, "vlnplot plot of TGFB1 and CXCL9 in tumor cells - remove zero values and log transformed experssion.pdf"),
  height = 10,
  width = 8
)

write.csv(
  dat_tgf_tumor_tidy2,
  row.names = FALSE,
  file  = file.path(analysis_out_dir, 'expression dat for TGFb.csv')
)


## T PR > PD_SD ---------
# p_vln_t_1 <- 
#   all_data2[, all_data2$group != "pre" & all_data2$annotation == "T"] %>%
#   VlnPlot(
#     features = c("CXCR3", "CXCR6"),
#     group.by = "group"
#   )  

dat_cxcr3_t <- 
  all_data2 %>%
  FetchData(vars = c("CXCR3", "CXCR6", "CCR5"), layer = "data") %>%
  cbind(all_data2[[c('annotation', 'celltype', 'group')]])


dat_cxcr3_t_tidy <-
  dat_cxcr3_t %>%
  dplyr::filter(group != "pre") %>% 
  dplyr::filter(annotation %in% "T") %>%
  pivot_longer(
    all_of(c("CXCR3", "CXCR6", "CCR5")),
    names_to = "gene",
    values_to = "expr"
  )  

dat_cxcr3_t_tidy2 <-
  dat_cxcr3_t_tidy %>%
  mutate(celltype = "T") %>% 
  rbind(dat_cxcr3_t_tidy) %>%
  
  # 20250630 adjust group
  mutate(
    group = fct_recode(
      group,
      Responder = "PR",
      'Non-responder' = 'PD_SD'
    ) %>%
      fct_relevel(
        c('Responder', 'Non-responder')
      )
  )

dat_cxcr3_t_tidy2$celltype <- factor(dat_cxcr3_t_tidy2$celltype, levels = c("T", t_cells))

p_cxcr3_t_1 <- 
  dat_cxcr3_t_tidy2 %>%
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(gene ~ celltype) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  # scale_y_continuous(transform = "log10") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) 


p_cxcr3_t_2 <- 
  dat_cxcr3_t_tidy2 %>%
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(gene ~ celltype) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  scale_y_continuous(transform = "log1p") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) 

p_cxcr3_t_3 <-
  dat_cxcr3_t_tidy2 %>%
  
  dplyr::filter(expr != 0) %>%   
  
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(gene ~ celltype) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  # scale_y_continuous(transform = "log1p") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) 

p_cxcr3_t_4 <-
  dat_cxcr3_t_tidy2 %>%
  
  dplyr::filter(expr != 0) %>%   
  
  ggplot(aes(x = group, y = expr)) +
  geom_violin(aes(fill = group), trim = FALSE) + 
  geom_boxplot(width = 0.15, fill = NA, outliers = FALSE) +
  facet_grid(gene ~ celltype) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  scale_y_continuous(transform = "log10") + # , limits = c(1e-3, 200)) +
  ggpubr::geom_pwc(
    tip.length = 0
  ) 

# export
ggsave(
  p_cxcr3_t_1,
  filename = file.path(analysis_out_dir, "vlnplot plot of CXCR3 CXCR6 and CCR5 in T cells.pdf"),
  height = 8,
  width = 25
)

ggsave(
  p_cxcr3_t_2,
  filename = file.path(analysis_out_dir, "vlnplot plot of CXCR3 CXCR6 and CCR5 in T cells - log transformed exression.pdf"),
  height = 8,
  width = 25
)

ggsave(
  p_cxcr3_t_3,
  filename = file.path(analysis_out_dir, "vlnplot plot of CXCR3 CXCR6 and CCR5 in T cells - remove zero values.pdf"),
  height = 8,
  width = 25
)

ggsave(
  p_cxcr3_t_4,
  filename = file.path(analysis_out_dir, "vlnplot plot of CXCR3 CXCR6 and CCR5 in T cells - remove zero values and log transformed experssion.pdf"),
  height = 8,
  width = 25
)


write.csv(
  dat_cxcr3_t_tidy2,
  row.names = FALSE,
  file  = file.path(analysis_out_dir, 'expression dat for CSCL_R.csv')
)

# 3. T cell score -------------
adata_t = sc$read('output/data/adata_t.h5ad')

expr_t <- t(as.matrix(adata_t$raw$X))
colnames(expr_t) <- adata_t$obs$cell_id_sample
rownames(expr_t) <- adata_t$raw$var_names$tolist()

meta_t <- adata_t$obs
rownames(meta_t) <- meta_t$cell_id_sample


exhaust_gene = c(
  'CTLA4',
  'PDCD1',
  'CXCL13',
  'ENTPD1',
  'LAG3',
  'LAYN',
  'TIGIT',
  'BATF',
  'HAVCR2',
  'TNFRSF9',
  'TOX'
)

active_gene = c(
  'KLRG1',
  'KLRD1',
  'GZMA',
  'GZMB',
  'GZMH',
  'GNLY',
  'ID2',
  'CCR2',
  'IFNG',
  'TBX21',
  'ENPTD1',
  'NKG7',
  'PRF1',
  'CCL5',
  'GPR183',
  'FGFBP2',
  'CX3CR1',
  'CD244',
  'SEMA4A',
  'GZMM',
  'FAS',
  'FASLG',
  'CD44',
  'CD69',
  'CD38',
  'NKG7',
  'KLRB1',
  'KLRD1',
  'KLRF1',
  'KLRG1',
  'KLRK1',
  'FCGR3A',
  'CX3CR1',
  'CD300A',
  'FGFBP2',
  'ID2',
  'ID3',
  'PRDM1',
  'RUNX3',
  'TBX21',
  'ZEB2',
  'BATF',
  'IRF4',
  'NR4A1',
  'NR4A2',
  'NR4A3',
  'PBX3',
  'ZNF683',
  'HOPX',
  'FOS',
  'FOSB',
  'JUN',
  'JUNB',
  'JUND',
  'STAT1',
  'STAT2',
  'STAT5A',
  'STAT6',
  'STAT4',
  'EOMES'
)

memory_gene <- c(
  'CCR7', 'CD28', 'CD44', 'CD69', 'CLDND1', 
  'CX3CR1', 'DKK3', 'DUSP2', 'EOMES', '
  GATA3', 'GPR183', 'GZMK', 'HSPA6', 'IL18RAP',
  'IL7R', 'ITGA1', 'KLF2', 'KLRG1,
  MT1E', 'MT1F', 'MYADM', 'NR4A1', 'PRDM1', 
  'RPS19', 'RPS26', 'TBX21,TGFBR2', 'TOB1', 
  'XCL1', 'XCL2', 'ZEB2'
)

# create seurat
not_pre_idx <- meta_t$group != 'pre'
all_data_t <-
  CreateSeuratObject(
    count = expr_t[, not_pre_idx],
    meta.data = meta_t[not_pre_idx, ]
  )

all_data_t@assays$RNA$data <- all_data_t@assays$RNA$counts

all_markers_t <-
  FindAllMarkers(
    all_data_t,
    group.by = "group"
  )


# exhaust gene: less expression in PR
all_markers_t_ex <-
  all_markers_t %>% 
  dplyr::filter(cluster == "PR") %>%
  dplyr::filter(gene %in% exhaust_gene) %>% 
  dplyr::filter(avg_log2FC < 0)

# active_gene: more expression in PR
all_markers_t_ac <-
  all_markers_t %>% 
  dplyr::filter(cluster == "PR") %>%
  dplyr::filter(gene %in% active_gene) %>% 
  dplyr::filter(avg_log2FC > 0)

# to python format
all_markers_t_ex$gene %>% str_c(collapse = "','") %>% {paste0("exhaust_gene2 = ['", ., "']")} %>% cat("\n")
all_markers_t_ac$gene %>% str_c(collapse = "','") %>% {paste0("active_gene2 = ['", ., "']")} %>% cat("\n")


# python socere --------------
geneset_ls = {
  'exhaust': exhaust_gene, 
  'active': active_gene,
  'exhaust2': exhaust_gene2, 
  'active2': active_gene2
}
for name, geneset in geneset_ls.items():
  sc.tl.score_genes(
    adata_t,
    gene_list=geneset,
    score_name=name
  )

# 4. T cell score - active and exhaust -------------
meta_t <- read.csv('report_all_samples/to_r_t_metadata.csv') 

p_t_geneset_ls <-
  c('exhaust', 'active', 'exhaust2', 'active2') %>%
  setNames(., .) %>%
  map(function(x){
    meta_t %>%
      dplyr::filter(group != "pre") %>%
      ggplot(aes(x = group, y = .data[[x]])) +
      geom_violin(aes(color = group)) +
      geom_boxplot(aes(color = group), width = 0.05, outliers = FALSE, show.legend = FALSE) +
      ggpubr::geom_pwc(tip.length = 0) +
      ggsci::scale_color_igv() +
      theme_classic() +
      ggtitle(x)
  })

# export 
pdf(
  file = file.path(str_glue("violin plot of geneset score - add pvalue.pdf")),
  height = 4,
  width  = 6 
)
p_t_geneset_ls %>%  walk(print)
dev.off()




#5.  percentage of memory gene  ----------
all_data_t_raw <- readRDS("report_all_samples/output/data/adata_t_raw.rds")
all_data_t_raw$nCount_RNA <- all_data_t_raw$total_counts

all_data_t_raw$percent.memory <- 
  all_data_t_raw %>%
  PercentageFeatureSet(features  = memory_gene)


p_vln_percent.memory_1 <- 
  all_data_t_raw[[]] %>% 
  dplyr::filter(group != "pre") %>% 
  dplyr::filter(percent.memory != 0) %>% 
  ggplot(aes(x = group, y = percent.memory, fill = group)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outliers = FALSE, show.legend = FALSE) +
  ggpubr::geom_pwc(tip.length = 0) +
  theme_classic() +
  ggtitle("remove zero values")


p_vln_percent.memory_2 <- 
  all_data_t_raw[[]] %>% 
  dplyr::filter(group != "pre") %>% 
  # dplyr::filter(percent.memory != 0) %>% 
  ggplot(aes(x = group, y = percent.memory, fill = group)) +
  geom_violin() +
  geom_boxplot(width = 0.1, outliers = FALSE, show.legend = FALSE) + 
  ggpubr::geom_pwc(tip.length = 0) +
  theme_classic() +
  ggtitle("all values")


# export 
pdf(
  file = file.path(analysis_out_dir, str_glue("violin plot of percent memory geneset.pdf")),
  height = 4,
  width  = 6 
)
p_vln_percent.memory_1 %>% print()
p_vln_percent.memory_2 %>% print()
dev.off()

# 6. barplot of cell props ---------

dat_props_group <- 
  all_data2[[]] %>%
  group_by(group, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) 


p_col_group <-
  dat_props_group %>%
  ggplot(aes(x = group, y = freq, fill = annotation)) +
  geom_col() +
  scale_fill_manual(values = cell_colors_used) +
  theme_classic() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )
  
dat_props_sample <- 
  all_data2[[]] %>%
  group_by(sample, annotation) %>%
  summarise(count = n()) %>%
  mutate(freq = count / sum(count)) 


p_col_sample <-
  dat_props_sample %>%
  ggplot(aes(x = sample, y = freq, fill = annotation)) +
  geom_col() +
  scale_fill_manual(values = cell_colors_used) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )


# export 
ggsave(
  p_col_group,
  filename = file.path(analysis_out_dir, "barplot of group.pdf"),
  height = 6, 
  width = 4
)

ggsave(
  p_col_sample,
  filename = file.path(analysis_out_dir, "barplot of sample.pdf"),
  height = 6, 
  width = 8
)

write.csv(
  dat_props_group,
  row.names = FALSE,
  file.path(analysis_out_dir, "cell props table - group.csv")
)

write.csv(
  dat_props_sample,
  row.names = FALSE,
  file.path(analysis_out_dir, "cell props table - sample.csv")
)


# 7. epxort cell annotation table -----------
dir.create("cell annotation table")

cell_annotation_ls <- 
  all_data[[]] %>% 
  dplyr::select(
    fov,
    cell_id	,
    group = celltype
  ) %>% 
  group_by(fov) %>%
  nest() %>%
  deframe() %>%
  map(as.data.frame)
  

cell_annotation_ls %>%
  iwalk(function(dat, name){
    write.csv(
      dat,
      row.names = FALSE,
      file.path('cell annotation table', str_glue("{name}.csv"))
    )
  })


# 8. tgfb expression ----------------

## plot boxplot and vlnplot -------
all_data2 %>%
  DotPlot(
    features = c("TGFB1", "TGFB1I1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3"),
    group.by = 'annotation'
  )

dat_tgfb <-
  all_data2 %>%
  FetchData(vars = c("TGFB1", "TGFB1I1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3")) %>%
  cbind(all_data2[[c('annotation', 'sample', 'group')]])

dat_tgfb_full <-
  dat_tgfb %>%
  pivot_longer(
    starts_with("TGF"),
    names_to = "TGFb",
    values_to = "expr"
  )

p_tgfb_ls_1 <-
  dat_tgfb_full %>%
  dplyr::filter(group != "pre") %>% 
  group_by(TGFb) %>%
  nest() %>%
  deframe() %>% # .[1] %>% 
  imap(function(dat, name){ # browser()
    dat %>%
      dplyr::filter(expr != 0) %>%
      dplyr::filter(annotation != "Epithelial") %>%
      ggplot(aes(x = annotation, y = expr, fill = group)) +
      # geom_violin() +
      geom_boxplot(outliers = FALSE) +
      geom_point(position = position_jitterdodge()) +
      ggpubr::geom_pwc(tip.length = 0) +
      # facet_wrap(~ TGFb, scales = "free_y") +
      # scale_y_continuous(transform = "log10")+
      theme_bw() +
      RotatedAxis() +
      labs(
        x = NULL, 
        y = 'log10(Normalized count)',
        title = name
      )
  })
  
p_tgfb_ls_2 <-
  dat_tgfb_full %>%
  dplyr::filter(group != "pre") %>% 
  group_by(TGFb) %>%
  nest() %>%
  deframe() %>% # .[1] %>% 
  imap(function(dat, name){
    dat %>%
      # dplyr::filter(expr != 0) %>%
      dplyr::filter(annotation != "Epithelial") %>%
      ggplot(aes(x = annotation, y = expr, color = group, fill = group)) +
      geom_violin() +
      geom_point(position = position_jitterdodge(), size = 0.5) +
      # geom_boxplot(outliers = FALSE) +
      # facet_wrap(~ TGFb, scales = "free_y") +
      # scale_y_continuous(transform = "log1p")+
      theme_bw() +
      RotatedAxis() +
      labs(
        x = NULL, 
        y = 'Normalized count',
        title = name
      )
  })
  

pdf(
  file = file.path(analysis_out_dir, "boxplot of TGFB gene in all major cells - remove zero values and log10 transformed experssion.pdf"),
  height = 6,
  width  = 12
)
p_tgfb_ls_1 %>% walk(print)
dev.off()
  

pdf(
  file = file.path(analysis_out_dir, "vlnplot of TGFB gene in all major cells.pdf"),
  height = 6,
  width  = 12
)
p_tgfb_ls_2 %>% walk(print)
dev.off()
  

##  plot heatmap  ------

p_ht_tgfb_ls <-
  dat_tgfb_full %>%
  
  dplyr::filter(group != "pre") %>%
  group_by(annotation, sample,   group, TGFb) %>%
  summarise(expr = mean(expr)) %>%
  
  group_by(TGFb) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  
  imap(function(dat, name){
    dat2 <- 
      dat %>%
      pivot_wider(
        names_from = "annotation",
        values_from = "expr"
      ) %>% 
      column_to_rownames('sample')
    
    ht <- 
      ComplexHeatmap::Heatmap(
      dat2[-1],
      name = "Expression",
      row_split = dat2$group,
      clustering_method_columns = "ward.D2",
      column_title = name,
      left_annotation = ComplexHeatmap::rowAnnotation(
        df = dat2[1], 
        col = list(group = c(PR = "#543005", PD_SD = "#003C30"))
      )
    )
    ht
  })
  
  


p_ht_tgfb_ls2 <-
  dat_tgfb_full %>%
  
  dplyr::filter(group != "pre") %>%
  group_by(annotation, sample,   group, TGFb) %>%
  summarise(expr = mean(expr)) %>%
  
  group_by(annotation) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  
  imap(function(dat, name){
    dat2 <- 
      dat %>%
      pivot_wider(
        names_from = "TGFb",
        values_from = "expr"
      ) %>% 
      column_to_rownames('sample')
    
    dat2[-1] <- scale(dat2[-1])
    val_rng_max <- range(dat2[-1]) %>% abs %>% max
    color_used <- RColorBrewer::brewer.pal(9, 'RdBu')
    ht <- 
      ComplexHeatmap::Heatmap(
        dat2[-1],
        name = "Expression",
        col = circlize::colorRamp2(
          colors = color_used,
          breaks = seq(-val_rng_max, val_rng_max, length.out = length(color_used))
        ),
        row_split = dat2$group,
        clustering_method_columns = "ward.D2",
        column_title = str_glue("{name} - scaled by column"),
        left_annotation = ComplexHeatmap::rowAnnotation(
          df = dat2[1], 
          col = list(group = c(PR = "#543005", PD_SD = "#003C30"))
        )
      )
    ht
  })
  
pdf(
  file = file.path(analysis_out_dir, "Heatmap of TGFB gene - in each TGFb gene.pdf"),
  height = 5,
  width  = 10
)
p_ht_tgfb_ls %>% walk(print)
dev.off()

pdf(
  file = file.path(analysis_out_dir, "Heatmap of TGFB gene - in each major celltype.pdf"),
  height = 4,
  width  = 8
)
p_ht_tgfb_ls2 %>% walk(print)
dev.off()


## plot heatmap - averaged all sample -----------
pdf(
  file = file.path(analysis_out_dir, "Heatmap of TGFB gene - average expr in each group and annotation.pdf"),
  height = 10,
  width  = 6
)
local({
  dat2 <-
    dat_tgfb_full %>%
  
    dplyr::filter(group != "pre") %>%
    group_by(annotation, group, TGFb) %>%
    summarise(expr = mean(expr)) %>%
    pivot_wider(
      names_from = "TGFb",
      values_from = "expr"
    ) 
  
  dat2[-(1:2)] <- scale(dat2[-(1:2)])
  val_rng_max <- range(dat2[-(1:2)]) %>% abs %>% max
  color_used <- RColorBrewer::brewer.pal(9, 'RdBu')
  
  ht <- 
    ComplexHeatmap::Heatmap(
      dat2[-(1:2)],
      name = "Expression",
      col = circlize::colorRamp2(
        colors = color_used,
        breaks = seq(-val_rng_max, val_rng_max, length.out = length(color_used))
      ),
      row_split = dat2$annotation,
      clustering_method_columns = "ward.D2",
      clustering_method_rows    = "ward.D2",
      row_title = " ",
      # column_title = str_glue(" "),
      # height = unit(8, "cm"),
      heatmap_height = unit(15, "cm"),
      heatmap_width = unit(8, "cm"),
      left_annotation = ComplexHeatmap::rowAnnotation(
        df = dat2[1:2], 
        col = list(
          group = c(PR = "#B15928", PD_SD = "#1F78B4"),
          annotation = cell_colors_used
        )
      )
    )
  ComplexHeatmap::draw(ht)
})
dev.off()


