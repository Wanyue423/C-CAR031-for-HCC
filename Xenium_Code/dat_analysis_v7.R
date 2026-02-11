library(tidyverse)
library(Seurat)

all_data2 <- readRDS('all_data2.rds')
all_data2$project <- all_data2$sample

# all_data2_cd45 <- readRDS('all_data2_cd45.rds')
# all_data2_t <- all_data2_cd45[, all_data2_cd45$annotation == "T"]


out_dir_v7 <- "20251113 dat analysis"
dir.create(out_dir_v7)




# pathway of immune signature ----------
pathway_ls <-  
  readxl::read_excel("13046_2018_1002_MOESM1_ESM.xlsx", skip = 1) %>% 
  map(~ .x[!is.na(.x)]) %>% 
  map(~ .x[.x %in% rownames(all_data2)]) %>% 
  discard(function(x) length(x) < 1) 


p_pathway_ls <- 
  pathway_ls %>%
  imap(function(x, name){ # browser()
    p_dot <- 
      DotPlot(
        all_data2,
        cols = c("navy", "white", "firebrick3"),
        group.by = "annotation",
        split.by = "group",
        features = x,
      )
    
    p_dot$data$annotation <- p_dot$data$id %>% str_remove("_.*$")
    p_dot$data$group <- p_dot$data$id %>% str_remove("^.*_")
    
    p <- p_dot$data %>%
      dplyr::filter(annotation %in% c("Macrophage", "Neutrophil")) %>% 
      dplyr::filter(group != "pre") %>%
      group_by(features.plot) %>%
      mutate(avg.exp.scaled = scale(avg.exp)) %>%  
      
      ggplot(aes(x = features.plot, y = group, size = pct.exp, color = avg.exp.scaled)) + 
      geom_point() +
      facet_wrap(~ annotation, ncol = 1) +
      theme_bw() +
      scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu")) ) +
      labs(
        title = name
      )
    
    # return
    list(
      p = p,
      len = length(x)
    )
  })




# 
# p_pathway_ls %>% 
#   iwalk(function(res, name){
#     p <- res$p
#     len <- res$len
#     
#     h <- 10 + len %/% 
#     
#   })


# export 
ggsave(
  p_pathway_ls$CCR$p + RotatedAxis(),
  filename = file.path(out_dir_v7, "1. mean expression in CCR.pdf"),
  height = 8,
  width  = 35
)

ggsave(
  p_pathway_ls$`Check-point`$p + RotatedAxis(),
  filename = file.path(out_dir_v7, "1. mean expression in checkpoint.pdf"),
  height = 8,
  width  = 12
)

ggsave(
  p_pathway_ls$`Inflammation-promoting`$p + RotatedAxis(),
  filename = file.path(out_dir_v7, "1. mean expression in Infla.pdf"),
  height = 8,
  width  = 12
)

ggsave(
  p_pathway_ls$Parainflammation$p + RotatedAxis(),
  filename = file.path(out_dir_v7, "1. mean expression in para-Infla.pdf"),
  height = 8,
  width  = 12
)

ggsave(
  p_pathway_ls$`APC co inhibition`$p + RotatedAxis(),
  filename = file.path(out_dir_v7, "1. mean expression in apc co-inhib.pdf"),
  height = 8,
  width  = 12
)

ggsave(
  p_pathway_ls$`APC co stimulation`$p + RotatedAxis(),
  filename = file.path(out_dir_v7, "1. mean expression in apc co-stimu.pdf"),
  height = 8,
  width  = 12
)




## selected gene 
p_dot <-
  DotPlot(
    all_data2,
    cols = c("navy", "white", "firebrick3"),
    group.by = "annotation",
    split.by = "group",
    features = c(
      'SERPINE1', 'THBS1', 'PDGFA',
      'VISTA', 'SIRPA', 'LILRB1', 'IL10',
      'SOCS1', 'ARG1',
      'CD74'
    )
  )

p_dot$data$annotation <- p_dot$data$id %>% str_remove("_.*$")
p_dot$data$group <- p_dot$data$id %>% str_remove("^.*_")


p_dot_selected  <-
  p_dot$data %>%
  dplyr::filter(annotation %in% c("Macrophage", "Neutrophil")) %>%
  dplyr::filter(group != "pre") %>%
  group_by(features.plot) %>%
  mutate(avg.exp.scaled = scale(avg.exp)) %>%

  ggplot(aes(x = features.plot, y = group, size = pct.exp, color = avg.exp.scaled)) +
  geom_point() +
  facet_wrap(~ annotation) +
  theme_bw() +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu")) ) +
  labs(
    title = ""
  )

# export 
ggsave(
  p_dot_selected + RotatedAxis(),
  filename = file.path(out_dir_v7, "3. mean expression in selected gene.pdf"),
  height = 8,
  width  = 12
)



# B cell signature ----------------
go_bp <- msigdbr::msigdbr(
  db_species = "HS",
  species = "human",
  collection = 'C5',
  subcollection = 'BP'
)

go_bp_b_cell_ls <-
  go_bp %>%
  dplyr::filter(str_detect(gs_name, "B_CELL")) %>%
  group_by(gs_name) %>% 
  nest() %>%
  deframe() %>%
  map('gene_symbol')


# de 
all_markers_b_cell_ls <-
  "B/Plasma cell" %>%
  setNames(., .) %>%
  map(function(x){
    Idents(all_data2) <- "annotation"
    all_markers <- 
      FindMarkers(
        all_data2,
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


# gsea
# all_gene_df <- 
#   clusterProfiler::bitr(
#     rownames(all_data2),
#     fromType = "SYMBOL",
#     toType = "ENTREZID",
#     OrgDb = "org.Hs.eg.db"
#   ) %>%
#   distinct() 


all_log2fc_b_cell_ls <-
  all_markers_b_cell_ls %>%
  map(function(all_marker){
    log2fc <- all_marker$avg_log2FC
    names(log2fc) <- all_marker$gene # all_gene_df$ENTREZID[match(all_marker$gene, all_gene_df$SYMBOL)]
    log2fc <- sort(log2fc, decreasing = TRUE)
    log2fc[!is.na(names(log2fc))]
  })



all_gsea_b_cell_full_ls <- 
  all_log2fc_b_cell_ls %>%
  map(function(x){
    clusterProfiler::GSEA(
      geneList = x,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500,
      TERM2GENE = go_bp_b_cell_ls %>% enframe("term", "gene") %>% unnest(gene)
    )
  })



p_dot_gsea_b_cell <-
  all_gsea_b_cell_full_ls$`B/Plasma cell`@result %>% 
  ggplot(aes(x = NES, y = reorder(Description, NES), fill = NES)) + 
  geom_col() + 
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu"))) + 
  theme( 
    legend.position = "none" 
  ) + 
  labs( 
    y = NULL 
  ) 


# export 
all_gsea_b_cell_full_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v7, "2. GSEA enrichment analysis - GOBP - only B cell pathway.xlsx"),
    asTable = TRUE
  )


ggsave(
  p_dot_gsea_b_cell,
  filename = file.path(out_dir_v7, "2. GSEA enrichment analysis  - GOBP - only B cell pathway - dotplot.pdf"),
  height = 8,
  width  = 8
)




# T cell signature ----------------
go_bp <- msigdbr::msigdbr(
  db_species = "HS",
  species = "human",
  collection = 'C5',
  subcollection = 'BP'
)

go_bp_t_cell_ls <-
  go_bp %>%
  dplyr::filter(str_detect(gs_name, "T_HELPER")) %>%
  group_by(gs_name) %>% 
  nest() %>%
  deframe() %>%
  map('gene_symbol')


# de 
all_markers_t_cell_ls <-
  "T" %>%
  setNames(., .) %>%
  map(function(x){
    Idents(all_data2) <- "annotation"
    all_markers <- 
      FindMarkers(
        all_data2,
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


# gsea
# all_gene_df <- 
#   clusterProfiler::bitr(
#     rownames(all_data2),
#     fromType = "SYMBOL",
#     toType = "ENTREZID",
#     OrgDb = "org.Hs.eg.db"
#   ) %>%
#   distinct() 


all_log2fc_t_cell_ls <-
  all_markers_t_cell_ls %>%
  map(function(all_marker){
    log2fc <- all_marker$avg_log2FC
    names(log2fc) <- all_marker$gene # all_gene_df$ENTREZID[match(all_marker$gene, all_gene_df$SYMBOL)]
    log2fc <- sort(log2fc, decreasing = TRUE)
    log2fc[!is.na(names(log2fc))]
  })



all_gsea_t_cell_full_ls <- 
  all_log2fc_t_cell_ls %>%
  map(function(x){
    clusterProfiler::GSEA(
      geneList = x,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 1,
      maxGSSize = 500,
      TERM2GENE = go_bp_t_cell_ls %>% enframe("term", "gene") %>% unnest(gene)
    )
  })



p_dot_gsea_t_cell <-
  all_gsea_t_cell_full_ls$`T`@result %>% 
  ggplot(aes(x = NES, y = reorder(Description, NES), fill = NES)) + 
  geom_col() + 
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(Inf, "RdBu"))) + 
  theme( 
    legend.position = "none" 
  ) + 
  labs( 
    y = NULL 
  ) 


# export 
all_gsea_b_cell_full_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_v7, "4. GSEA enrichment analysis - GOBP - only B cell pathway.xlsx"),
    asTable = TRUE
  )


ggsave(
  p_dot_gsea_t_cell,
  filename = file.path(out_dir_v7, "4. GSEA enrichment analysis  - GOBP - only T cell pathway - dotplot.pdf"),
  height = 8,
  width  = 8
)













