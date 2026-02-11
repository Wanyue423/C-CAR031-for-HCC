rm(list = ls())

load(file = "Tumor_Step2.Rdata")
getwd()
gmt_h <- readLines("./GSEA/h.all.v2024.1.Hs.symbols.gmt")
gmt_c4<-readLines("./GSEA/c4.all.v2024.1.Hs.symbols.gmt")
gmt_c7<-readLines("./GSEA/c7.all.v2024.1.Hs.symbols.gmt")

gene_set1 <- strsplit(gmt_h[grepl("HALLMARK_INTERFERON_GAMMA_RESPONSE", gmt_h)], "\t")[[1]][-c(1,2)]
gene_set2 <- strsplit(gmt_h[grepl("HALLMARK_OXIDATIVE_PHOSPHORYLATION", gmt_h)], "\t")[[1]][-c(1,2)]
gene_set3 <- strsplit(gmt_h[grepl("HALLMARK_MYC_TARGETS_V1", gmt_h)], "\t")[[1]][-c(1,2)]
gene_set4 <- strsplit(gmt_h[grepl("HALLMARK_ANGIOGENESIS", gmt_h)], "\t")[[1]][-c(1,2)]

a = ""
b = ""
c = ""
gene_sig1=a
gene_sig2=b
gene_sig3=c

data_gene_set <- deg %>%
  dplyr::filter(symbol%in%c(gene_set1, gene_set2,gene_set3,gene_set4,gene_sig1,gene_sig2)) %>%
  dplyr::mutate(gene_set=case_when(
    symbol%in%gene_sig1~"gene_sig1",
    symbol%in%gene_sig2~"gene_sig2",
    symbol%in%gene_set1~"gene_set1",
    symbol%in%gene_set2~"gene_set2",
    symbol%in%gene_set3~"gene_set3",
    symbol%in%gene_set4~"gene_set4"))
data_other_genes<-deg%>%
  filter(!symbol%in%c(gene_set1,gene_set2,gene_set3,gene_set4,gene_sig1,gene_sig2))

down_genes<-deg%>%
  dplyr::filter(P.Value<0.05,logFC < -0.5)
up_genes<-deg%>%
  dplyr::filter(P.Value<0.05,logFC > 0.5)
top_down200_genes<-deg%>%
  filter(P.Value<0.01,logFC<0)%>%
  arrange(logFC)%>%
  head(200)
top_up200_genes<-deg%>%
  filter(P.Value<0.01,logFC>0)%>%
  arrange(desc(logFC))%>%
  head(200)

p=ggplot(deg,aes(x=logFC,y=-log10(P.Value)))+
  annotate("rect",xmin=min(up_genes$logFC),xmax=2.86,ymin=-log10(0.05),
           ymax=25.3,fill= "#FBB4AE")+
  annotate("rect",xmin=-2.90,xmax=max(down_genes$logFC),ymin=-log10(0.05),
           ymax=25.3,fill="#A6CEE3")+
  annotate("rect",xmin=min(top_up200_genes$logFC),xmax=2.86,ymin=-log10(0.01),
           ymax=25.3,fill="transparent",linetype="dotted",color="black",linewidth=0.6)+
  annotate("rect",xmin=-2.90,xmax=max(top_down200_genes$logFC),ymin=-log10(0.01),
           ymax=25.3,fill="transparent",linetype="dotted",color="black",linewidth=0.6)+

geom_vline(xintercept=0,lty=4,lwd=0.6,alpha=0.8)+
  geom_hline(yintercept=0,lty=4,lwd=0.6,alpha=0.8)+
  geom_hline(yintercept=-log10(0.05),lty=4,lwd=0.6,alpha=0.8)+
  geom_point(data=data_other_genes,shape=21,color="black",alpha=0.2,size=1.2,stroke=0.7)+
  geom_point(data=data_gene_set,aes(fill=gene_set),shape=21,color="black",size=1.8,stroke=0.8)+
  geom_label_repel(data = filter(data_gene_set, gene_set == "gene_sig1"),
                   aes(label = symbol), size = 6, col = "black", 
                   box.padding = unit(0.5, 'lines'), fontface = "bold", fill = "white",alpha = 0.8,max.overlaps = 10,nudge_y = 0,          # 微调标签在y轴方向的偏移
                   nudge_x = 0) +
  
  geom_label_repel(data = filter(data_gene_set, gene_set == "gene_sig2"),
                   aes(label = symbol), size = 6, col = "black", 
                   box.padding = unit(0.5, 'lines'), fontface = "bold", fill = "white",alpha = 0.8,max.overlaps = 10,nudge_y = -0.2,          # 微调标签在y轴方向的偏移
                   nudge_x = 0)+
  
geom_text_repel(
  data = filter(data_gene_set, gene_set == "gene_sig3"),
  aes(label = symbol),
  size = 3,
  color = "black",  
  fontface = "bold",
  max.overlaps = 10,
  nudge_y = 0,
  nudge_x = 0,
  box.padding = unit(0.5, "lines"),  
  segment.color = "grey50",          
  min.segment.length = 0.2          
)+
  scale_y_continuous(
    name = expression(-log[10]("P value")),
    limits = c(-12, 31),
    breaks = seq(0, 32, 6),
    expand = c(0, 0)
  )+
  scale_x_continuous(name = expression(log[2] * "FC"), limits = c(-3, 3), breaks = seq(-3,3,3), expand = c(0, 0.1)) +
  
  scale_fill_manual(values=c("gene_set1"="#62b448","gene_set2"= "#3e56a1","gene_set3"="#9f9ea3","gene_set4"="#9F815B","gene_sig1"="#07F1F9","gene_sig2"="#FF0000","gene_sig3"="transparent"))+
  
  theme(panel.border = element_rect(fill=NA,color="black", size = 0.8,linetype="solid"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size = 14,color="black"),
        axis.text = element_text(size = 28,color="black"),
        axis.title.x = element_text(size = 28,color="black"),
        axis.title.y = element_text(size = 28,color="black",vjust = 2),
        axis.line = element_blank()) +
  
  geom_linerange(data=filter(data_gene_set,gene_set=="gene_set1"),aes(x=logFC,ymin=-2.7,ymax=-0.3),
                 color="#62b448",size=13,linewidth=0.3)+
  geom_linerange(data=filter(data_gene_set,gene_set=="gene_set2"),aes(x=logFC,ymin=-5.7,ymax=-3.3),
                 color="#3e56a1",size=13,linewidth=0.3)+
  geom_linerange(data=filter(data_gene_set,gene_set=="gene_set3"),aes(x=logFC,ymin=-8.7,ymax=-6.3),
                 color="#9f9ea3",size=13,linewidth=0.3)+
  geom_linerange(data=filter(data_gene_set,gene_set=="gene_set4"),aes(x=logFC,ymin=-11.7,ymax=-9.3),
                 color="#9F815B",size=13,linewidth=0.3)+
  
  
  
  annotate("text",x=-1.57,y=28.2,label="Differentially-regulated",color="#1F78B4",size=11,lineheight=0.8,vjust=0,fontface = "bold")+
  annotate("text",x=1.54,y=28.2,label="Differentially-regulated",color= "#E31A1C",size=11,lineheight=0.8,vjust=0,fontface = "bold")+
  annotate("text",x=-2.69,y=25.8,label="(PR)",color="#1F78B4",size=10,lineheight=0.8,vjust=0,fontface = "plain")+
  annotate("text",x=2.39,y=25.9,label="(Non-PR)",color= "#E31A1C",size=10,lineheight=0.8,vjust=0,fontface = "plain")+
  annotate("text",x=-2.37,y=3.5,label="Top 200",color="black",size=11,fontface = "bold")+
  annotate("text",x=2.30,y=3.5,label="Top 200",color="black",size=11,fontface = "bold")+
  annotate("text",x=-3,y=-1.5,label="INTERFERON_GAMMA_RESPONSE",color="#62b448",size=9.5,hjust=0,fontface = "bold")+
  annotate("text",x=-3,y=-4.5,label="OXIDATIVE_PHOSPHORYLATION",color="#3e56a1",size=9.5,hjust=0,fontface = "bold")+
  annotate("text",x=-3,y=-7.5,label="MYC_TARGETS_V1",color="#9f9ea3",size=9.5,hjust=0,fontface = "bold")+
  annotate("text",x=-3,y=-10.5,label="ANGIOGENESIS",color="#9F815B",size=9.5,hjust=0,fontface = "bold")+
  
  annotate("text",x=3.0,y=-1.5,label="NES:-2.31 pval:5.9e-13",color="#62b448",size=9,lineheight=0.8,hjust=1,fontface = "bold")+
  annotate("text",x=3.0,y=-4.5,label="NES:2.28 pval:2.0e-11",color="#3e56a1",size=9,lineheight=0.8,hjust=1,fontface = "bold")+
  annotate("text",x=3.0,y=-7.5,label="NES:2.04 pval:1.6e-07",color="#9f9ea3",size=9,lineheight=0.8,hjust=1,fontface = "bold")+
  annotate("text",x=3.0,y=-10.5,label="NES:-1.80 pval:2.5e-03",color="#9F815B",size=9,lineheight=0.8,hjust=1,fontface = "bold")

p

