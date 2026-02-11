library(tidyverse)
library(Seurat)

all_data2 <- readRDS('all_data2.rds')
all_data2$project <- all_data2$sample

out_dir_revised <- "20250924"
dir.create(out_dir_revised)

cell_colors_used <- 
  all_data2$celltype %>%
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

# define function ------------
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

# Immune cell props ---------
all_data2$annotation %>% unique
# [1] "Tumor cell"                        "Hepatocyte"                        "Macrophage"                       
# [4] "Epithelial"                        "Endothelial"                       "Proliferating tumor cell"         
# [7] "T"                                 "Kupffer cell"                      "Cholangiocyte"                    
# [10] "Pericyte"                          "B/Plasma cell"                     "Hepatic stellate cell"            
# [13] "Liver sinusoidal endothelial cell" "Neutrophil"                        "Macrophage-CAF"  


cd45_cells <- 
  c(
    "Macrophage",
    "T",
    "Kupffer cell",
    "B/Plasma cell",
    "Neutrophil"
  )

p_immune_cell_pie <-
  all_data2[[]] %>%
  mutate(
    annotation2 = ifelse(
      annotation %in% cd45_cells,
      "Immune cells",
      as.character(annotation)
    )
  ) %>%
  dplyr::filter(
    annotation2 != "Macrophage-CAF",
    group != "pre"
  ) %>%
  group_by(group, annotation2) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(
    freq = count/sum(count),
    group = factor(group, levels = c("PD_SD", "PR")),
    # annotation2 = str_c(annotation2, scales::percent(freq, accuracy = 0.1), sep = " - ")
  ) %>%

  ggplot(aes(x = 1, y = freq, fill = annotation2)) +
  geom_col() + 
  # ggrepel::geom_text_repel(
  #   data = ~ mutate(., freq2 = cumsum(freq), freq3 = (c(0, freq2[-length(freq2)]) + freq2)/2),
  #   aes(y = freq3, label = signif(freq, 2))
  # ) +
  coord_polar(theta = "y") +
  facet_wrap(~ group) +
  theme_void() +
  theme(plot.margin = margin(1, 1, 1, 1, "mm"))  +
  scale_fill_manual(
    values = c(cell_colors_used, "Immune cells" = "#E7475E")
  )

# p_immune_cell_pie

ggsave(
  p_immune_cell_pie, 
  filename = file.path(out_dir_revised, "1 pie plot of immune cells props.pdf"),
  height = 5,
  width = 8
)


# umap plot of all immune cells props ----------

all_data2_cd45 <- 
  all_data2[,
    all_data2$annotation %in% cd45_cells &
      (!str_detect(all_data2$celltype, "like")) &
      (all_data2$group != "pre")
  ]

all_data2_cd45 <-
  all_data2_cd45 %>%
  # NormalizeData(scale.factor = median(.$nCount_RNA)) %>% # logNormalize在r/python中是相同的
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = 'total_counts') %>%
  RunPCA()

ElbowPlot(all_data2_cd45, ndims = 30)


dims <- 10


# library(future)
# plan(multisession, workers = 10)
# options(future.globals.maxSize = 10 * 1024 ^ 3) # 10Gb


all_data2_cd45 <- RunUMAP(all_data2_cd45, dims = 1:dims, return.model = F)


DimPlot(all_data2_cd45, group.by = "project", label = TRUE, reduction = "umap") + NoLegend()


all_data2_cd45 <- all_data2_cd45 %>%
  harmony::RunHarmony(
    dims.use = 1:dims,
    group.by.vars = "project",
    reduction.save = "harmony"
  )

all_data2_cd45 <- RunUMAP(
  all_data2_cd45,
  dims = 1:dims,
  reduction = "harmony",
  reduction.name = "harmonyumap",
  reduction.key  = "harmonyUMAP_", 
  return.model = TRUE
)


# save
saveRDS(all_data2_cd45, file = 'all_data2_cd45.rds')

DimPlot(all_data2_cd45, group.by = "annotation", label = TRUE, reduction = "umap") + NoLegend()
DimPlot(
  all_data2_cd45, 
  group.by = "annotation", 
  label = TRUE, 
  reduction = "harmonyumap"
) + 
  NoLegend() +
  scale_color_manual(values = col_plan_1)


df_freq_immune_cd45 <-
  all_data2_cd45[[]] %>%
  dplyr::filter(group != "pre") %>%
  group_by(group, annotation) %>%
  summarise(count = n()) %>%
  group_by(group) %>%
  mutate(freq = count/sum(count))
  
  
df_freq_immune_cd45_toadd <-
  all_data2_cd45[[c('group', 'annotation')]] %>%
  rownames_to_column('id') %>%
  dplyr::filter(group != "pre") %>%
  left_join(df_freq_immune_cd45, by = c("group", "annotation")) %>%
  mutate(annotation_tidy = paste0(annotation, " (", scales::percent(freq, accuracy = 0.1), ")")) %>%
  column_to_rownames('id') %>%
  .['annotation_tidy']
  
col_for_immune_subcluster <- c(
  'Macrophage (37.0%)' = "#ea7070",
  'Macrophage (37.9%)' = "#ea7070",

  'T (34.1%)'  = "#fdc4b6",
  'T (25.9%)'  = "#fdc4b6",

  'Kupffer cell (6.2%)'  = "#e59572",
  'Kupffer cell (15.9%)' = "#e59572",

  'B/Plasma cell (4.1%)' = "#2694ab",
  'B/Plasma cell (4.4%)' = "#2694ab",

  'Neutrophil (18.6%)'   = "#96ceb4",
  'Neutrophil (15.9%)'   = "#96ceb4"
)


p_umap_immune_subcluster <-
  all_data2_cd45[, all_data2_cd45$group != "pre"] %>%
  AddMetaData(df_freq_immune_cd45_toadd) %>%
  DimPlot(
    group.by = "annotation_tidy", 
    label = TRUE, 
    reduction = "harmonyumap",
    split.by = "group"
  ) + 
  NoLegend()  +
  scale_color_manual(values = col_for_immune_subcluster)

ggsave(
  add_minicoord(p_umap_immune_subcluster), 
  filename = file.path(out_dir_revised, "2 umap plot of immune subcluster cells.pdf"),
  height = 5,
  width = 10
)

# export cell freq by sample

all_data2[[]] %>%
  dplyr::filter(
    all_data2$annotation %in% cd45_cells,
    !str_detect(all_data2$celltype, "like")
  )  %>%
  # dplyr::filter(group != "pre") %>%
  group_by(sample, annotation) %>%
  summarise(count = n()) %>%
  group_by(sample) %>%
  mutate(freq = count/sum(count)) %>%
  dplyr::select(-count) %>%
  pivot_wider(
    names_from = "annotation",
    values_from = "freq"
  ) %>%
  
  openxlsx::write.xlsx(
    file = file.path(out_dir_revised, "2 cell props of immune subcluster cells.xlsx"),
    asTable = TRUE
  )

# de analysis of immune subcluster between group ------------
Idents(all_data2_cd45) <- 'celltype'

all_markers_cd45_ls <-
  all_data2_cd45$celltype %>%
  unique() %>%
  setNames(., .) %>%
  map(function(x){
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

 
all_gene_cd45_ls <-
  all_markers_cd45_ls %>%
  map(function(all_marker){
      up_gene <- 
        dplyr::filter(
          all_marker,
          avg_log2FC > 0.25,
          p_val_adj < 0.05
        ) %>%
        dplyr::pull('gene')
      
      down_gene <- 
        dplyr::filter(
          all_marker,
          avg_log2FC < -0.25,
          p_val_adj < 0.05
        ) %>%
        dplyr::pull('gene')
      
      list(
        PR = up_gene,
        PD_SD = down_gene
      )
  })
  
# some t cell have no significant genes
# using gsva
all_gene_df <- 
  clusterProfiler::bitr(
    rownames(all_data2_cd45),
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  ) %>%
  distinct() %>%
  dplyr::filter()

all_log2fc_cd45_ls <-
  all_markers_cd45_ls %>%
  map(function(all_marker){
    log2fc <- all_marker$avg_log2FC
    names(log2fc) <- all_gene_df$ENTREZID[match(all_marker$gene, all_gene_df$SYMBOL)]
    log2fc <- sort(log2fc, decreasing = TRUE)
    log2fc[!is.na(names(log2fc))]
  })



all_gsea_cd45_ls <- 
  all_log2fc_cd45_ls %>%
  map(function(x){
    clusterProfiler::gseGO(
      gene = x,
      ont = "BP",
      OrgDb = "org.Hs.eg.db",
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
  })

p_dot_gsea_cd_45 <-
  all_gsea_cd45_ls %>%
  imap(function(enrich, name){
    tryCatch(
      enrichplot::dotplot(
        enrich,
        showCategory = 20,
        label_format = 50,
        title = name
      ),
      error = function(e) return(NULL)
    )
    
  })

p_gsea_plot_cd45 <-
  all_gsea_cd45_ls %>%
  imap(function(enrich, name){
    p_blank <- ggplot()
    paths <- all_gsea_cd45_ls$Macrophage@result$Description[1:5]
    1:5 %>%
      map(function(x){
        tryCatch(
          enrichplot::gseaplot2(enrich, geneSetID = x, title = paths[x]),
          error = function(e) return(p_blank)
        )
      })
  })

# export
pdf(
  file = file.path(out_dir_revised, "3 dot plot of immune cells.pdf"),
  height = 6,
  width = 8
)
p_dot_gsea_cd_45 %>% walk(print)
dev.off()


p_gsea_plot_cd45 %>%
  iwalk(function(p, name){
    name <- str_replace(name, "/", "_")
    pdf(
      file = file.path(out_dir_revised, str_glue("3 gsea plot of immune cells - {name}.pdf")),
      height = 6,
      width = 8
    )
    p %>% walk(print)
    dev.off()
  })

# export 
all_gsea_cd45_ls %>%
  map('result') %>%
  openxlsx::write.xlsx(
    file = file.path(out_dir_revised, "3 gsea enrichment table.xlsx"),
    asTable = TRUE
  )


