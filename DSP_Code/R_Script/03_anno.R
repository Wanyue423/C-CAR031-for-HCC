rm(list = ls())  
load(file = 'Tumor_Step2.Rdata')

gene_up = deg[deg$change == 'up','ENTREZID'] 
gene_down = deg[deg$change == 'down','ENTREZID'] 
gene_diff = c(gene_up,gene_down)

ego_diff <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "ALL",
                     readable = TRUE)

ego_diff <- as.data.frame(ego_diff)

write.csv(ego_diff, "ego_diff.csv", row.names = FALSE)
