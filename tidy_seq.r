tidy_seq <- function (x) {
  if(is.null(x) == F){
    stop("x should not be empty.")
}
 # retrieving to official name 
 # BioMart database 
  mart <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
 # extracting with matching hgnc symbol
 gene_official_name <- biomaRt::getBM(values = x,
                                filters = "hgnc_symbol",
                                mart = mart,
                                attributes = c("external_gene_name", "hgnc_symbol",
                                               "entrezgene_id","ensembl_gene_id", 
                                               "uniprotswissprot","ensembl_transcript_id" ))
 

 
 transcriptCoordsByGene.GRangesList_unlist <-GenomicFeatures::transcriptsBy( TxDb.Hsapiens.UCSC.hg38.knownGene,
                                                            by = "exon") %>%
   plyranges::as_granges() %>% 
   plyranges::mutate( tx_name = gsub("\\..*","",tx_name)) %>% 
   plyranges::filter (tx_name %in% unique(gene_official_name$ensembl_transcript_id))

  trs_data <- keepStandardChromosomes(transcriptCoordsByGene.GRangesList_unlist,species="Homo_sapiens",
                                   pruning.mode="coarse")
  
  return(trs_data)
}
