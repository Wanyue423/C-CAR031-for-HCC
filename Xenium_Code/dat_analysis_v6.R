library(tidyverse)
library(Seurat)

all_data2 <- readRDS('all_data2.rds')
all_data2$project <- all_data2$sample

# all_data2_cd45 <- readRDS('all_data2_cd45.rds')
# all_data2_t <- all_data2_cd45[, all_data2_cd45$annotation == "T"]


out_dir_v6 <- "20251106 dat analysis"
dir.create(out_dir_v6)


# WNT module score ------------
pathway_wnt_manual_ls <-
  readxl::excel_sheets("KEGG_WNT_pathway_filter_by_expr.xlsx") %>%
  map(~{
    readxl::read_excel("KEGG_WNT_pathway_filter_by_expr.xlsx", sheet = .x)
  }) 

pathway_wnt_manual_ls2 <- 
  pathway_wnt_manual_ls[[1]] %>% 
  group_by(pathway) %>%
  nest() %>%
  deframe() %>%
  map('SYMBOL') %>%
  map(c, pathway_wnt_manual_ls[[2]]$SYMBOL)


pathway_wnt_manual_ls2[['other']] <- c(
  pathway_wnt_manual_ls2$PCP,
  pathway_wnt_manual_ls2$`wnt/ca`
)

# all_data2_sub <- all_data2[rowMeans(all_data2[['RNA']]$data) > 0.01, ]

all_data2_score <- all_data2 %>% AddModuleScore(features = pathway_wnt_manual_ls2)
# colnames(all_data2_score@meta.data)[
#   str_detect(colnames(all_data2_score@meta.data), "Cluster")
# ] <- names(pathway_wnt_manual_ls2)
  


df_wnt_manual <-
  all_data2_score %>%
  FetchData(
    vars = c(
      str_c("Cluster", seq_along(pathway_wnt_manual_ls2)),
      'celltype', 'annotation',
      'group',
      'project'
    )
  )


df_wnt_manual_tidy <- 
  df_wnt_manual %>%
  group_by(project, group, celltype) %>%
  dplyr::select(-annotation) %>%
  summarise(across(everything(), mean)) %>%
  dplyr::filter(
    group != "pre",
    celltype == "Tumor cell"
  ) 


df_heat_wnt_manual <- 
  df_wnt_manual_tidy %>%
  
  pivot_longer(
    starts_with("Cluster"),
    names_to = "pathway",
    values_to = "val"
  ) %>%
  
  group_by(pathway) %>%
  mutate(val = scales::rescale(val, to = c(0,1))) %>%
  
  group_by(group, celltype, pathway) %>%
  summarise(val = mean(val)) 


p_heat_wnt_manual <- 
  df_heat_wnt_manual %>% 
  
  mutate(pathway = fct_recode(pathway, !!!structure(
    names = names(pathway_wnt_manual_ls2),
    str_c("Cluster", seq_along(pathway_wnt_manual_ls2))
  ))) %>% 
  
  ggplot(aes(x = pathway, y = group, fill = val)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "MinMax normalized score"
  )

# export
pdf(
  file = file.path(out_dir_v6, 'WNT pathway heatmap.pdf'),
  height = 3.3,
  width  = 5
)
p_heat_wnt_manual %>% print()
dev.off() 

df_heat_wnt_manual %>% 
  
  mutate(pathway = fct_recode(pathway, !!!structure(
    names = names(pathway_wnt_manual_ls2),
    str_c("Cluster", seq_along(pathway_wnt_manual_ls2))
  ))) %>%
  write.csv(
    file = file.path(out_dir_v6, 'WNT pathway heatmap value.csv'),
    row.names = FALSE
  )



# TGFb module score ------------

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





# TGFb score
all_data2_score2 <- all_data2 %>% AddModuleScore(features = h_geneset_ls['tgf_beta_signaling'])

p_vln_tgf <-
  all_data2_score2 %>%
  VlnPlot(
   features = 'Cluster1',
   group.by = 'annotation',
   split.by = 'group',
   pt.size = 0
  ) +
  facet_wrap(~ split)


df_heat_tgf <-
  all_data2_score2[[c('group', 'annotation', 'Cluster1')]] %>%
  # group_by(group) %>%
  # mutate(val = scales::rescale(val, to = c(0,1))) %>%
  # 
  group_by(group, annotation) %>%
  summarise(Cluster1 = mean(Cluster1)) 


p_heat_tgf  <-
  df_heat_tgf %>% 
  
  group_by(annotation) %>%
  mutate(Cluster1.scaled = scale(Cluster1)) %>%
  
  ggplot(aes(x = annotation, y = group, fill = Cluster1.scaled)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("#010107", "#BD3583", "#FCFFA4"))(50)) +
  theme_minimal() +
  theme(
    legend.position = "top",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "TGFb pathway score"
  ) +
  guides(
    fill = guide_colorbar(theme = theme(
      legend.key.width  = unit(10, "lines")
    ))
  )




# export
ggsave(
  p_heat_tgf,
  filename = file.path(out_dir_v6, 'TFGb pathway heatmap.pdf'),
  height = 3.3,
  width  = 8
)


ggsave(
  p_vln_tgf,
  filename = file.path(out_dir_v6, 'TFGb pathway vlnplot.pdf'),
  height = 6,
  width  = 12
)

df_heat_tgf %>% 
  dplyr::rename(TGFb = Cluster1) %>%
  write.csv(
    file = file.path(out_dir_v6, 'TGFb pathway heatmap value.csv'),
    row.names = FALSE
  )





p_dot_tgf <- 
  DotPlot(
    all_data2_score2,
    group.by = "annotation",
    features = "Cluster1",
    split.by = "group",
     cols = c("navy", "white",  "firebrick3")
  )


# p_dot_tgf$data %>% head


# all_data2_score2[[c('annotation', 'Cluster1')]]  %>%
#   group_by(annotation) %>%
#   summarise(mean(expm1(Cluster1)) %>% log1p())






p_dot_cell_grp <- 
  p_dot_tgf$data %>%
  mutate(
    cell = str_remove(id, "_.*$"),
    group  = str_remove(id, "^.*_")
  ) %>%
  ggplot(aes(x = group, y = cell, color = avg.exp, size = pct.exp)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds") ) +
  theme_bw() +
  labs(
    size = "pct.exp.gt.0"
  )

p_dot_cell_grp2 <- 
  p_dot_tgf$data %>%
  mutate(avg.exp.scaled = scale(avg.exp)) %>%
  mutate(
    cell = str_remove(id, "_.*$"),
    group  = str_remove(id, "^.*_")
  ) %>%
  ggplot(aes(x = group, y = cell, color = avg.exp.scaled, size = pct.exp)) +
  geom_point(shape = 16) +
  scale_color_gradientn(colours = RColorBrewer::brewer.pal(Inf, "Reds") ) +
  theme_bw() +
  labs(
    size = "pct.exp.gt.0"
  )




ggsave(
  p_dot_cell_grp,
  filename = file.path(out_dir_v6, 'TFGb pathway dotplot.pdf'),
  height = 6,
  width  = 6
)
ggsave(
  p_dot_cell_grp2,
  filename = file.path(out_dir_v6, 'TFGb pathway dotplot2.pdf'),
  height = 6,
  width  = 6
)







