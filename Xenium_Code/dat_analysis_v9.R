library(tidyverse)
library(Seurat)

all_data2_cd45 <- readRDS('all_data2_cd45.rds')

all_data2 <- readRDS('all_data2.rds')
all_data2$project <- all_data2$sample

# all_data2_cd45 <- readRDS('all_data2_cd45.rds')
# all_data2_t <- all_data2_cd45[, all_data2_cd45$annotation == "T"]

all_data_epi <- readRDS('all_data_epi.rds')

out_dir_v9 <- "20251207 dat analysis"
dir.create(out_dir_v9)


# 1. GPC3+ cells proportion ----------------
pre_grp <- c(
  PFS_long = "CJU-pre",
  PFS_long = "HYLO-pre",
  PFS_long = "XXLI-pre",
  PFS_long = "LBCH-pre",
  
  PFS_short = "XWRO-pre",
  PFS_short = "LZQI-pre",
  PFS_short = "LXZH-pre"
)


 
p_col_pfs <-
  all_data_epi[[]] %>%  
  dplyr::filter(sample %in% pre_grp) %>% 
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>% 
  group_by(PFS_group, AUC_GPC3_Group) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = AUC_GPC3_Group )) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  labs(
    y = NULL
  )

# export 
pdf(
  file.path(out_dir_v9, "1 GPC3 cells - cell prop of PFS group.pdf"), 
  height = 8, 
  width  = 5
)
p_col_pfs %>% cowplot::plot_grid() %>% print()
dev.off()

## 1.2 by patient ------------

p_box_pfs <-
  all_data_epi[[]] %>%  
  dplyr::filter(sample %in% pre_grp) %>% 
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>% 
  
  group_by(sample, PFS_group, AUC_GPC3_Group) %>%
  summarise(count = n()) %>%
  
  group_by(sample, PFS_group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = PFS_group )) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitter(width = 0.35)) +
  ggpubr::geom_pwc(tip.length = 0) +
  facet_wrap(~ AUC_GPC3_Group) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  labs(
    y = NULL
  )


# export 
pdf(
  file.path(out_dir_v9, "1 GPC3 cells - cell prop of PFS group - by sample - boxplot.pdf"), 
  height = 5, 
  width  = 8
)
p_box_pfs %>% print()
dev.off()


# 2. nn cell props of PFS long and PFS shor ------------
res_nn_ls <- readRDS("20251022 epi nn analysis/res_nn_ls.rds")

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
p_nn_props_col_pfs_only_gpc3 <-
  res_nn_name_df_gpc3_full %>%
  
  dplyr::filter(
    project %in% pre_grp,
    AUC_GPC3_Group == 'GPC3'
  ) %>% 
  mutate(PFS_group = fct_recode(project, !!! pre_grp)) %>% 
  
  group_by(PFS_group, NN_cells) %>%
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = PFS_group, y = freq)) +
  geom_col(aes(fill = NN_cells)) +
  
  scale_fill_manual(values = col_plan_2) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = NULL,
    y = NULL
  )


ggsave(
  p_nn_props_col_pfs_only_gpc3,
  filename = file.path(out_dir_v9, "2 nn cells - bar plot of gpc3 - group by PFS group.pdf"), 
  height = 8, 
  width  = 5
)


## 2.2 by patient ------------
# plot
p_nn_props_box_pfs_only_gpc3 <-
  res_nn_name_df_gpc3_full %>%
  
  dplyr::filter(
    project %in% pre_grp,
    AUC_GPC3_Group == 'GPC3'
  ) %>% 
  mutate(PFS_group = fct_recode(project, !!! pre_grp)) %>% 
  
  group_by(project, PFS_group, NN_cells) %>%
  summarise(count = n()) %>%
  
  group_by(project, PFS_group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = PFS_group)) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitter(width = 0.35)) +
  ggpubr::geom_pwc(tip.length = 0, bracket.nudge.y = 0.1) +
  facet_wrap(~ NN_cells, scales = "free_y") +
  scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]) ) +
  theme_classic() +
  theme(
    strip.background = element_blank()
  ) +
  labs(
    y = NULL
  )


ggsave(
  p_nn_props_box_pfs_only_gpc3,
  filename = file.path(out_dir_v9, "2 nn cells - boxplot of gpc3 - group by PFS group - by sample.pdf"), 
  height = 12, 
  width  = 15
)




# 3. cell props of PFS group -------------
p_col_pfs_group <-
  all_data2[[]] %>%
  dplyr::filter(sample %in% pre_grp) %>%
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>%
  group_by(PFS_group, annotation) %>%
  
  summarise(count = n()) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = annotation)) +
  geom_col() +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = col_plan_2) +
  labs(
    y = NULL
  )


ggsave(
  p_col_pfs_group,
  filename = file.path(out_dir_v9, "3 cell prop of PFS group.pdf"), 
  height = 8, 
  width  = 5
)

## 3.2 by patient ------------ 
p_box_pfs_group <-
  all_data2[[]] %>%
  dplyr::filter(sample %in% pre_grp) %>%
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>%
  
  group_by(sample, PFS_group, annotation) %>%
  summarise(count = n()) %>%
  
  group_by(sample, PFS_group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = PFS_group )) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitter(width = 0.35)) +
  ggpubr::geom_pwc(tip.length = 0, bracket.nudge.y = 0.1) +
  facet_wrap(~ annotation, scales = "free_y") +
  theme_classic() +
  scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]) ) +
  theme(
    strip.background = element_blank()
  ) +
  labs(
    y = NULL
  )


ggsave(
  p_box_pfs_group,
  filename = file.path(out_dir_v9, "3 cell prop of PFS group - by sample - boxplot.pdf"), 
  height = 12, 
  width  = 15
)



## 3.3 by patient - T / total cells ----------
# 20251209

t_subtype <- c(
  'Tregs',
  'Hepatocyte-like T cell',
  'Gamma delta T',
  'Macrophage-like T cell',
  'CD4T',
  'CD8Tem',
  'Fibroblast-like T cell',
  'Cycling T cell'
)

p_box_pfs_t_group <-
  all_data2[[]] %>%
  dplyr::filter(sample %in% pre_grp) %>%
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>%
  
  group_by(sample, PFS_group, celltype) %>%
  summarise(count = n()) %>%
  
  group_by(sample, PFS_group) %>%
  mutate(freq = count/sum(count)) %>%
  
  dplyr::filter(celltype %in% t_subtype) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = PFS_group )) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitter(width = 0.35)) +
  ggpubr::geom_pwc(tip.length = 0, bracket.nudge.y = 0.1) +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_classic() +
  scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]) ) +
  theme(
    strip.background = element_blank()
  ) +
  labs(
    y = NULL
  )


ggsave(
  p_box_pfs_t_group,
  filename = file.path(out_dir_v9, "3 cell prop of PFS group - by sample - boxplot - t cells among all cells.pdf"), 
  height = 9, 
  width  = 12
)




## 3.4 by patient - T / immune cells ----------
# 20251209
p_box_pfs_t_cd45_group <-
  all_data2[[]] %>%
  dplyr::filter(
    annotation %in% c(
       "Macrophage",
       "T",
       "Kupffer cell",
       "B/Plasma cell",
       "Neutrophil"
     )
  ) %>% 
  dplyr::filter(sample %in% pre_grp) %>%
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>%
  
  group_by(sample, PFS_group, celltype) %>%
  summarise(count = n()) %>%
  
  group_by(sample, PFS_group) %>%
  mutate(freq = count/sum(count)) %>%
  
  dplyr::filter(celltype %in% t_subtype) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = PFS_group )) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitter(width = 0.35)) +
  ggpubr::geom_pwc(tip.length = 0, bracket.nudge.y = 0.1) +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_classic() +
  scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]) ) +
  theme(
    strip.background = element_blank()
  ) +
  labs(
    y = NULL
  )


ggsave(
  p_box_pfs_t_cd45_group,
  filename = file.path(out_dir_v9, "3 cell prop of PFS group - by sample - boxplot - t cells among CD45 cells.pdf"), 
  height = 9, 
  width  = 12
)


## 3.5 by patient - T / immune cells ----------
# 20251209
p_box_pfs_t_t_group <-
  all_data2[[]] %>%
  dplyr::filter(annotation == "T") %>% 
  dplyr::filter(sample %in% pre_grp) %>%
  mutate(PFS_group = fct_recode(sample, !!! pre_grp)) %>%
  
  group_by(sample, PFS_group, celltype) %>%
  summarise(count = n()) %>%
  
  group_by(sample, PFS_group) %>%
  mutate(freq = count/sum(count)) %>%
  
  dplyr::filter(celltype %in% t_subtype) %>%
  
  ggplot(aes(x = PFS_group, y = freq, fill = PFS_group )) +
  geom_boxplot(outliers = F) +
  geom_point(position = position_jitter(width = 0.35)) +
  ggpubr::geom_pwc(tip.length = 0, bracket.nudge.y = 0.1) +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_classic() +
  scale_y_continuous(limits = function(x) c(x[1], 1.2 * x[2]) ) +
  theme(
    strip.background = element_blank()
  ) +
  labs(
    y = NULL
  )


ggsave(
  p_box_pfs_t_t_group,
  filename = file.path(out_dir_v9, "3 cell prop of PFS group - by sample - boxplot - t cells among t cells.pdf"), 
  height = 9, 
  width  = 12
)
