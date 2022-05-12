wTO_FIMO <- function(dataset, adjustpvalue = 0.05 , 
                     wTOscore = 0.5, TFlist ,
                     Node1= NULL, Node2 =NULL,
                     upstream = 5e3, downstream= 3e3 ) {
  # checking for missing input
  # dataset should be represented  in hgnc symbol
  if(is.data.frame(dataset) == F){
    stop("dataset must be a data.frame.")
  }
  if(is.numeric(adjustpvalue) == F){
    stop("adjustpvalue must be numeric.")
  }
  if(wTOscore <= 0){
    stop("adjustpvalue must be greater than 0.")
  }
  if(is.numeric(wTOscore) == F){
    stop("adjustpvalue must be numeric.")
  }
  if(wTOscore <= 0){
    stop("wTOscore must be greater than 0.")
  }
  if(is.data.frame(TFlist) == 0){
    stop("TFlist must be a data.frame.")
  }
  if(is.character(TFlist) == F){
    stop("TFlist must be a data.frame.")
  }
   if(is.null(Node1) == F & is.null(Node2) == F & is.null(TFlist) == F){
    stop("TF list or Node 1 and Node 2 should not be empty")
  }

 # filtering data   
 dataset_use<- dataset %>% dplyr::filter(pval.adj<adjustpvalue, abs(wTO) >= wTOscore) 
 unique_genelist<- unique(c(dataset_use$Node.1,dataset_use$Node.2))
 
clean_sequence<- tidy_seq(unique_genelist)
promter_seq <- retrieve_seq(clean_sequence,upstream,downstream)
thres<- (abs(upstream) +abs(downstream))/1000 # every n*1kb, n will be divided
result <- FIMO(promter_seq, TFlist, thres)
 
rm(thres)
rm(dataset_use)
rm(unique_genelist)
return(result)
} 
