options("repos"="https://mirrors.ustc.edu.cn/CRAN/")
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

cran_packages <- c('tidyr',
                   'tibble',
                   'dplyr',
                   'stringr',
                   'ggplot2',
                   'ggpubr',
                   'factoextra',
                   'FactoMineR',
                   'devtools',
                   'cowplot',
                   'patchwork',
                   'basetheme',
                   'paletteer') 
Biocductor_packages <- c('GEOquery',
                         'hgu133plus2.db',
                         'ggnewscale',
                         "limma",
                         "impute",
                         "GSEABase",
                         "GSVA",
                         "clusterProfiler",
                         "org.Hs.eg.db",
                         "preprocessCore",
                         "enrichplot",
                         "ggplotify")

for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    require(pkg,character.only=T) 
  }
}


for (pkg in c(Biocductor_packages,cran_packages)){
  require(pkg,character.only=T) 
}

#If there is no prompt, it means success. If there is a warning that a certain package is missing, check it using library.
#If an error occurs, go back and reinstall. If you haven't installed a certain package but are prompted that it is missing, this is normal due to complex dependencies. Just install whatever is missing.

if(!require(AnnoProbe))devtools::install_local("./AnnoProbe-master.zip",upgrade = F)
library(AnnoProbe)
library(openxlsx)
library(sva)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(tidyverse)
library(limma)
library(readxl)
