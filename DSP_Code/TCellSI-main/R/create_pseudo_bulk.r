#' @title TCSS_Calculate
#' @description You can use this function to calculate scores for T cell states.
#'
#' @param object Sample expression
#' @param reference Reference inputs
#' @param nbin How many gene bins
#' @param ctrl How many genes to select from each bins
#' @param seed Random number seed set
#' @return features.scores.df
#' @import ggplot2
#' @import dplyr
#' @importFrom stats na.omit rnorm
#' @export
#' @examples
#'
#' # Example usage:
#' \donttest{
#' sample_expression <- TCellSI::exampleSample
#' TCSS_Calculate(sample_expression)
#' }
rm(list = ls())
load("CD3_combat_count.Rdata")
library(TCellSI)
ResultScores <- TCellSI::TCSS_Calculate(exp)

create_pseudo_bulk <- function(annotation_data, expression_data, cluster_col, cell_id_col, n_clusters = 18, factor = 5, sampling_rate = 0.6) {
  # 计算伪散装数据列表
  pseudo_bulk_list <- lapply(1:n_clusters, function(i) {
    pseudo_bulk_df <- as.data.frame(lapply(1:round(table(annotation_data[[cluster_col]])[i] / factor), function(j) {
      sp <- sample(annotation_data[which(annotation_data[[cluster_col]] == names(table(annotation_data[[cluster_col]]))[i]), ][[cell_id_col]], round(table(annotation_data[[cluster_col]])[i] / factor) * sampling_rate)
      pseudo_bulk <- as.data.frame(apply(expression_data[, which(colnames(expression_data) %in% sp)], 1, mean))
      return(pseudo_bulk)
    }))
    return(pseudo_bulk_df)
  })
  pseudo_bulk <- as.data.frame(matrix(nrow = nrow(expression_data), ncol = 0))
  for (i in 1:n_clusters) {
    colnames(pseudo_bulk_list[[i]]) <- rep(paste0(names(table(annotation_data[[cluster_col]]))[i], "_bulk"), round(table(annotation_data[[cluster_col]])[i] / factor))
    pseudo_bulk <- cbind(pseudo_bulk, pseudo_bulk_list[[i]])
  }
  return(pseudo_bulk)
}

Ttype= ResultScores
p = identical(rownames(pd),colnames(Ttype));p
if (!p) {
  s = intersect(rownames(pd), colnames(Ttype))
  Ttype = Ttype[, s, drop = FALSE]  # 保持 exp 的数据框结构
  pd = pd[s, , drop = FALSE]     # 确保 pd 仍然是数据框
}

score=as.data.frame(t(Ttype))
dat1 = t(Ttype) %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  mutate(group = Group)
write.xlsx(dat1, file = "Ttype.xlsx")
library(FactoMineR)
library(factoextra)
dat.pca <- PCA(score, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = Group, # color by groups
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups"
)
g = names(tail(sort(apply(exp,1,sd)),1000))
n = exp[g,]
library(pheatmap)
annotation_col = data.frame(row.names = colnames(n),
                            Group = Group)
library(limma)
design=model.matrix(~Group)
fit=lmFit(Ttype,design)
fit=eBayes(fit)
score_final=topTable(fit,coef=2,number = Inf)


