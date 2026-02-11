rm(list = ls()) 

setwd("/Users/whipple1/Desktop/Bioinfomatics/DSP_DEG")
pd <- read.xlsx("Group_PD_PR.xlsx", rowNames = TRUE) 
exp <- read.xlsx("Expression_q_norm.counts.xlsx", rowNames = TRUE) 
pd <- na.omit(pd)
k = str_detect(pd$PD_PR,"Tumor_PR|Tumor_PD|Tumor_Pre");table(k)
pd = pd[k, , drop = FALSE]
Group=pd$PD_PR
Group = factor(Group,levels = c("Tumor_PR","Tumor_PD","Tumor_Pre"))
p = identical(rownames(pd),colnames(exp));p
if (!p) {
  s = intersect(rownames(pd), colnames(exp))
  exp = exp[, s, drop = FALSE]  
  pd = pd[s, , drop = FALSE]     
}


pre.pca <- PCA(t(exp),graph = FALSE)
fviz_pca_ind(pre.pca,
             geom= c("point"),
             col.ind = pd$slide,
             addEllipses = TRUE,
             legend.title="Group")

combat_exp <- ComBat_seq(counts = as.matrix(exp), 
                         batch = pd$slide,
                         group=Group)

combat.pca <- PCA(t(combat_exp),graph = FALSE)
fviz_pca_ind(combat.pca,col.ind = pd$slide,
             geom = c("point"),
             addEllipses = TRUE,
             legend.title="Group" )
exp=log2(combat_exp+ 1)
save(pd,Group,exp,file = "Tumor_combat_count.Rdata")

