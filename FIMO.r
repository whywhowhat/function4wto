
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


FIMO <- function(promter_seq , TFlist, thres  ){
if(is.null(promter_seq) == F){
    stop("promter_seq should not be empty.")
}if(is.character(TFlist) == F){
    stop("thres must be character")
}if(is.numeric(thres) == F){
    stop("thres must be numeric.")
  }
  
Fimo_result_datalist <- list()
  for (i in 1:length(TFlist)) {  #HOCOMOCO$`Transcription factor`
   tryCatch({
     getting_moitf<- MotifDb::MotifDb %>%
       # Query the database for the motif using it's gene name
       MotifDb::query(TFlist[i],andStrings= c("hsapiens","hocomoco")) %>% #HOCOMOCO$`Transcription factor`
       # Convert from motifdb format to universalmotif format
       universalmotif::convert_motifs() %>%
       # The result is a list, to simplify the object, return it as a single universalmotif
       .[[1]]
     #running fimo
     Fimo_result <- runFimo(promter_seq, getting_moitf,
                            silent = T,
                            thresh = 1e-4/round(thres,1)
                            )

     # print(i) # optional to see progress
     Fimo_result_dataframe <- Repitools::annoGR2DF(Fimo_result)   # grange to dataframe
     Fimo_result_dataframe$i <- i # keeping track of the iteration
     Fimo_result_datalist[[i]] <- Fimo_result_dataframe # adding to list
   },error=function(e){cat("Node",
                           TFlist[i],
                           "not found :",conditionMessage(e), "\n")})
} %>% suppressMessages()


Combind_data <- dplyr::bind_rows(Fimo_result_datalist)
Combind_data <- Combind_data %>%
  dplyr::mutate(set_of_motif_seq = str_c(Combind_data$motif_id,
                                         Combind_data$matched_sequence,
                                  sep = "_and_"))


return(Combind_data)
}
