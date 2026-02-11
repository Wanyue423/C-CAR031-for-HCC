rm(list = ls()) 
load(file = "Tumor_combat_count.Rdata")

design <- model.matrix(~0 + Group) 
colnames(design) <- levels(Group)   
rownames_identical <- identical(colnames(exp), rownames(pd))

contrast.matrix <- makeContrasts(
  PD_vs_PR = Tumor_PD - Tumor_PR,
  levels = design
)

fit <- lmFit(exp, design)       
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

PD_vs_PR=topTable(fit2, coef = "PD_vs_PR", number = Inf)

logFC_t=0.5
P.Value_t = 0.05
k1 = (PD_vs_PR$P.Value < P.Value_t)&(PD_vs_PR$logFC < -logFC_t)
k2 = (PD_vs_PR$P.Value < P.Value_t)&(PD_vs_PR$logFC > logFC_t)
deg <- mutate(PD_vs_PR,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
table(deg$change)

deg <- mutate(deg, symbol = rownames(deg))
s1e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
deg <- inner_join(deg,s1e,by=c("symbol"="SYMBOL"))

deg <- deg[!grepl("RPL|RPS", deg$symbol), ]
save(deg,logFC_t,P.Value_t,file = "Tumor_Step2.Rdata")

