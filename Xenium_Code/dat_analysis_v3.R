library(tidyverse)
library(Seurat)

all_data2_cd45 <- readRDS('all_data2_cd45.rds')

out_dir_v3 <- "20251009 pathway score of non-T cell"
dir.create(out_dir_v3)

# pathway of immune signature ----------
pathway_ls <- 
  readxl::read_excel("13046_2018_1002_MOESM1_ESM.xlsx", skip = 1) %>%
  map(~ .x[!is.na(.x)]) %>%
  map(~ .x[.x %in% rownames(all_data2_cd45)]) %>%
  discard(function(x) length(x) < 1)


all_data2_cd45_score <-
  all_data2_cd45 %>%
  AddModuleScore(pathway_ls)


df_pathway <- 
  all_data2_cd45_score[[]] %>%
  dplyr::select(starts_with("Cluster")) %>%
  setNames(names(pathway_ls))


df_pathway_tidy <- 
  all_data2_cd45_score[[c('group', 'sample', 'annotation')]]  %>%
  cbind(df_pathway) %>%
  as.data.frame() %>%
  dplyr::filter(annotation != "T", group != "pre") 
  
  
p_pathway_group_boxplot_ls <-
  df_pathway_tidy %>% 
  pivot_longer(
    -c(sample, group, annotation),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  group_by(sample, group, annotation, pathway) %>%
  summarise(score = mean(score)) %>%
  
  group_by(pathway) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  
  imap(function(dat, name){
    dat %>%
      ggplot(aes(x = group, y = score, fill = group)) +
      geom_boxplot(outliers = FALSE) +
      geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
      facet_wrap(~ annotation) +
      theme_bw() +
      labs(
        title = name
      )
  })
  

# export
pdf(
  file.path(out_dir_v3, "pathway score of non-T cell.pdf"),
  height = 6,
  width  = 8
)
p_pathway_group_boxplot_ls %>% walk(print)
dev.off()

  
# pathway of myeloid signature 1 ----------
mye_pathway_ls <- 
  readxl::read_excel("41467_2023_36296_MOESM5_ESM.xlsx", skip = 2, sheet = 'Figure 5') %>%
  map(~ .x[!is.na(.x)]) %>%
  map(~ .x[.x %in% rownames(all_data2_cd45)]) %>%
  discard(function(x) length(x) < 1)


all_data2_cd45_mye1 <-
  all_data2_cd45 %>%
  AddModuleScore(mye_pathway_ls)


df_pathway_mye1 <- 
  all_data2_cd45_mye1[[]] %>%
  dplyr::select(starts_with("Cluster")) %>%
  setNames(names(mye_pathway_ls))


df_pathway_tidy_mye1 <- 
  all_data2_cd45_mye1[[c('group', 'sample', 'annotation')]]  %>%
  cbind(df_pathway_mye1) %>%
  as.data.frame() %>%
  dplyr::filter(annotation != "T", group != "pre") 


p_mye1_pathway_group_boxplot_ls <-
  df_pathway_tidy_mye1 %>% 
  pivot_longer(
    -c(sample, group, annotation),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  group_by(sample, group, annotation, pathway) %>%
  summarise(score = mean(score)) %>%
  
  group_by(pathway) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  
  imap(function(dat, name){
    dat %>%
      ggplot(aes(x = group, y = score, fill = group)) +
      geom_boxplot(outliers = FALSE) +
      geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
      facet_wrap(~ annotation) +
      theme_bw() +
      labs(
        title = name
      )
  })


# export
pdf(
  file.path(out_dir_v3, "pathway score of non-T cell - myeloid signature 1.pdf"),
  height = 6,
  width  = 8
)
p_mye1_pathway_group_boxplot_ls %>% walk(print)
dev.off()

df_pathway_tidy_mye1 %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v3, "pathway score of non-T cell - myeloid signature 1.xlsx"),
    asTable = TRUE
  )

# pathway of myeloid signature 2 ----------
mye_pathway_ls2 <- 
  readxl::read_excel("pahtway of myeloid.xlsx", skip = 0) %>%
  map(~ .x[!is.na(.x)]) %>%
  map(str_to_upper) %>%
  map(~ .x[.x %in% rownames(all_data2_cd45)]) %>%
  discard(function(x) length(x) < 1)
mye_pathway_ls2 <- mye_pathway_ls2[!str_detect(names(mye_pathway_ls2), "AT2")]

all_data2_cd45_mye2 <-
  all_data2_cd45 %>%
  AddModuleScore(mye_pathway_ls2)


df_pathway_mye2 <- 
  all_data2_cd45_mye2[[]] %>%
  dplyr::select(starts_with("Cluster")) %>%
  setNames(names(mye_pathway_ls2))


df_pathway_tidy_mye2 <- 
  all_data2_cd45_mye2[[c('group', 'sample', 'annotation')]]  %>%
  cbind(df_pathway_mye2) %>%
  as.data.frame() %>%
  dplyr::filter(annotation != "T", group != "pre") 


p_mye2_pathway_group_boxplot_ls <-
  df_pathway_tidy_mye2 %>% 
  pivot_longer(
    -c(sample, group, annotation),
    names_to = "pathway",
    values_to = "score"
  ) %>%
  group_by(sample, group, annotation, pathway) %>%
  summarise(score = mean(score)) %>%
  
  group_by(pathway) %>%
  nest() %>%
  deframe() %>% # .[1] %>%
  
  imap(function(dat, name){
    dat %>%
      ggplot(aes(x = group, y = score, fill = group)) +
      geom_boxplot(outliers = FALSE) +
      geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.25)) +
      facet_wrap(~ annotation) +
      theme_bw() +
      labs(
        title = name
      )
  })


# export
pdf(
  file.path(out_dir_v3, "pathway score of non-T cell - myeloid signature 2.pdf"),
  height = 6,
  width  = 8
)
p_mye2_pathway_group_boxplot_ls %>% walk(print)
dev.off()

df_pathway_tidy_mye2 %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v3, "pathway score of non-T cell - myeloid signature 2.xlsx"),
    asTable = TRUE
  )


