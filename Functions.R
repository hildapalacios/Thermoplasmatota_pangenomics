get_expanded_classification <- function(in_DF,colname){
  input_classification <- in_DF[,grep(colname,colnames(in_DF))] %>% dplyr::pull(colname) %>% unique() # generates a non redundant input vector for get_ids
  results <- taxize::classification(input_classification, db = "ncbi", rows = 1, messages = FALSE) # save results of non interactive get_ids
  results_unlisted <- results %>% unlist() # unlist results
  all_ranks <- results_unlisted[grep("rank", names(results_unlisted))] %>% unname() %>% unique() # complete repertoire of taxonomic ranks
  
  out_DF <- matrix(nrow = length(results), ncol = length(all_ranks) + 1) %>% as.data.frame(.) # output data frame many rows as organisms, many columns as ranks
  colnames(out_DF) <- c(colname, paste(colname, c(all_ranks), sep = "_")) # naming of taxonomic rank columns
  
  # for loop to match taxonomic classifications
  index = 0 # iterative variable
  for (e in input_classification) {
    index = index + 1 # update iterative variable
    clasif <- results[[e]] # retrieve classification found for organism e
    out_DF[index ,which(colnames(out_DF) == colname)] <- e # write query name in first column
    if(any(is.na(clasif))){next}
    # for each species, iterate over each global taxonomic range
    for (r in all_ranks) {
      r_classif <- clasif$name[grep(r,clasif$rank)] # given lineage for rank r
      r_outDF_index <- grep(r, colnames(out_DF)) # define column index of rank r
      if(length(r_classif)==0){out_DF[index,r_outDF_index] <- ""}else{out_DF[index,r_outDF_index] <- r_classif[1]} # conditionally fill output df if a given rank was assigned for a given query
    }
  }
  return(out_DF)
}
