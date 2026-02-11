library(tidyverse)
library(Seurat)

all_data2_cd45 <- readRDS('all_data2_cd45.rds')

out_dir_v4 <- "20251016"
dir.create(out_dir_v4)


# 1. immune cell / CD45 cells --------
all_data2_cd45[[]] %>%
  group_by(project, group, annotation) %>%
  summarise(count = n()) %>%
  group_by(project, group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = annotation, y = freq, color = group, fill = group)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75)) +
  ggpubr::geom_pwc(tip.length = 0, method = "wilcox.test") +
  theme_bw()


all_data2_cd45[[]] %>%
  group_by(project, group, celltype) %>%
  summarise(count = n()) %>%
  group_by(project, group) %>%
  mutate(freq = count/sum(count)) %>%
  
  ggplot(aes(x = celltype, y = freq, color = group, fill = group)) +
  geom_boxplot(alpha = 0.5, outliers = FALSE) +
  geom_point(position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75)) +
  ggpubr::geom_pwc(tip.length = 0, method = "wilcox.test") +
  theme_bw()



#2   Ro/e -------
anno_column <- "annotation"
grp_column  <- "project"

all_data_meta_roe <-
  all_data2_cd45[[]] %>%
  mutate(
    {{anno_column}} := fct_drop(.data[[anno_column]]),
    {{grp_column}}  := fct_drop(.data[[grp_column]])
  )



# Cell Preference
library(ComplexHeatmap)
# 生成细胞数目分布频数统计表observe.data
observe.data <- unclass(table(all_data_meta_roe[[anno_column]], all_data_meta_roe[[grp_column]]))
stopifnot(is.matrix(observe.data))

expected <-
  function(x){
    if(!is.matrix(x)){stop("Must be a matrix")}
    rtot <- margin.table(x, 1)
    ctot <- margin.table(x, 2)
    tot  <- margin.table(x)
    outer(rtot, ctot, "*")/tot
  }
expected.data <- expected(observe.data) # 生成期望频数表
plot.data <- observe.data/expected.data # 计算实际细胞数目与期望细胞数目的比值
# 绘图数据的列顺序
# col.order <- ggplot2:::`%||%`(
#   levels(all_data_meta_roe[[grp_column]]),
#   unique(all_data_meta_roe[[grp_column]])
# ) 

col.order <- all_data2_cd45[[c('group', "sample")]] %>% distinct() %>% arrange(group) %>% pull('sample') %>% as.character()

plot.data <- plot.data[, col.order] # 对绘图数据列进行排序
plot.data <- plot.data[order(apply(plot.data, 1, which.max), decreasing = T), ] # 对绘图数据行进行排序
# 如果想深色色块在主对角线使用which.min，在副对角线则使用which.max

# 热图颜色设定
roe_cols <- structure(
  c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
  names = c("1", "1.5", "2", ">3")
)

cell_fun = function(j, i, x, y, width, height, fill) { # browser()
  
  # print(c(j, i, x, y, width, height, fill))
  
  if(plot.data[i,j] < as.numeric(names(roe_cols)[1])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[1], fill = roe_cols[1]))
    grid::grid.text("\u00B1", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[2])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[2], fill = roe_cols[2]))
    grid::grid.text("+", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(roe_cols)[3])) {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[3], fill = roe_cols[3]))
    grid::grid.text("++", x = x, y = y)
  } else {
    grid::grid.rect(x = x, y = y, width = width, height = height, gp = grid::gpar(col = roe_cols[4], fill = roe_cols[4]))
    grid::grid.text("+++", x = x, y = y)
  }
  
}

# 生成热图
roi_hm <-
  ComplexHeatmap::Heatmap(
    plot.data, 
    cell_fun = cell_fun,
    cluster_rows = F, 
    cluster_columns = F,
    width = unit(8, "cm"), 
    height = unit(11, "cm"), 
    show_heatmap_legend = F
  )

# stat comparasion
chi2.T <- sum((observe.data - expected.data)^2/expected.data)
p_chisq <- pchisq(q = chi2.T, df = (nrow(observe.data)-1)*(ncol(observe.data)-1), lower.tail = FALSE)
print(p_chisq)

# 生成图例
lgd <-
  ComplexHeatmap::Legend(
    labels = names(roe_cols), 
    title = expression(R[o/e]), 
    legend_gp = grid::gpar(fill = roe_cols)
  )

# pdf(
#   file.path(roe_out_dir, "cell props boxplot - Roe plot.pdf"), 
#   width = 8, 
#   height = 6
# )
ComplexHeatmap::draw(roi_hm, annotation_legend_list = lgd)
# dev.off()






