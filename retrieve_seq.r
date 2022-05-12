retrieve_seq <- function (x,upstream ,downstream){
  if(is.null(x) == F){
    stop("x should not be empty.")
}if(is.numeric(upstream) == F){
    stop("upstream must be numeric.")
}if(is.numeric(downstream) == F){
    stop("downstream must be numeric.")
  }
 # retrieving promoter seqeunce from
 promoter_seqs <-  x %>% 
   plyranges::anchor_start() %>% # strech seqeunce to upstream
   plyranges::stretch(upstream) %>% 
   plyranges::anchor_end() %>% # strech seqeunce to  downstream
   plyranges::stretch(downstream) %>% 
   get_sequence(BSgenome.Hsapiens.UCSC.hg38)
 
 return(promoter_seqs)
}
