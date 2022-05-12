# function4wto
to enrich wTO package 






# installing package seperately 
install.packages("dplyr")
install.packages("magrittr")
install.packages("stringr")
BiocManager::install("GenomicFeatures",force = TRUE  )
BiocManager::install("GenomicRanges",force = TRUE  )
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene",force = TRUE ) 
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38",force = TRUE )
BiocManager::install("MotifDb",force = TRUE)
BiocManager::install("universalmotif",force = TRUE)
BiocManager::install("biomaRt",force = TRUE)
BiocManager::install("plyranges",force = TRUE)
BiocManager::install("Repitools")#


suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressPackageStartupMessages(library(MotifDb))
suppressPackageStartupMessages(library(universalmotif))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(Repitools))
suppressPackageStartupMessages(library(memes))
# should localy installed
library(memes)
