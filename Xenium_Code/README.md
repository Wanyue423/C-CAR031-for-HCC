# GPC3-specific dnTGFβRII–armoured CAR-T cells for hepatocellular carcinoma

### pathway score and Differential analysis  

Xenium-derived single-cell expression matrices and spatial coordinates were analyzed using Seurat (v5). 
Pathway activity was quantified using gene set–based scoring. TGFβ signaling activity was calculated using the Hallmark TGFβ gene set from MSigDB with AddModuleScore.
Differential gene expression was performed using FindMarkers in macrophages, CD4⁺ T cells, and B/plasma cells. Gene set enrichment analysis (GSEA) was conducted using clusterProfiler, with enrichment quantified by normalized enrichment scores. Hallmark and KEGG pathways were used for macrophages, and Gene Ontology Biological Process terms were used for CD4⁺ T cells and B/plasma cells.


### Cell neighborhood   


Tumour spatial neighbourhoods were analyzed using a k-nearest neighbour (kNN) approach. For each tumour cell, the 100 nearest neighbouring cells were identified based on Euclidean spatial distances, and neighbourhood cell-type proportions were calculated. The resulting neighbourhood matrix was subjected to PCA and graph-based clustering using Seurat, with the first five principal components retained.


### GPC3 cells identification    


GPC3-positive tumour cells were inferred by projecting a GPC3-associated gene signature, derived from paired DSP differential expression analysis, onto the Xenium dataset using AUCell. Mixture modeling of AUCell score distributions was used to define classification thresholds, and cells exceeding the threshold were classified as GPC3^high.


### other codes    

Other codes used to reproduce figures in this paper.



