#install.packages("pacman")
pacman::p_load(BiocManager, tidyverse)
#BiocManager::install(c("Rfastp"))
library(Rfastp)

# functions
createFilePairsList <- function(paths, pattern) {
  common_parts <- gsub(paste0(pattern,"\\.fastq\\.gz"), "", paths) %>% basename()
  unique_common_parts <- unique(common_parts)
  file_pairs_list <- list()
  for(common_part in unique_common_parts){
    matched_files <- grep(common_part, paths, value = TRUE)
    if(length(matched_files) == 2 )　{
      file_pairs_list[[common_part]] <- list(
        common_part = common_part,
        file1 = matched_files[1],
        file2 = matched_files[2])
    }
  }
  return(file_pairs_list)
}

extract_summary_parts <- function(youso) {
  # 'before_filtering' と 'after_filtering' 部分を抽出し、データフレームに変換
  before_filtering_df <- youso$summary$before_filtering %>% as.data.frame() %>% t() %>% data.frame() %>% rownames_to_column()
  after_filtering_df <- youso$summary$after_filtering %>% as.data.frame() %>% t() %>% data.frame() %>% rownames_to_column()
  names(before_filtering_df) = c("item", "before_filtering")
  names(after_filtering_df) = c("item", "after_filtering")
  return( inner_join(before_filtering_df, after_filtering_df, by="item"))
}

extract_whatever_filtering <- function(youso){
  total_reads <- youso$total_reads
  total_bases <- youso$total_bases
  q20_bases <- youso$q20_bases
  q30_bases <- youso$q30_bases
  quality_curves <- youso$quality_curves %>% as.data.frame()
  content_curves <- youso$content_curves %>% as.data.frame()
  kmer_count <- youso$kmer_count %>% as.data.frame() %>% t()
  return(list(main=data.frame(total_reads, total_bases, q20_bases, q30_bases), quality_curves = quality_curves, content_curves = content_curves, kmer_count = kmer_count))  
}


# Rfastpの実行
runRfastp <- function(result_dir_parsed, file_pair) {
  rfastp_result <- Rfastp::rfastp(read1 = file_pair$file1, 
                                  read2 = file_pair$file2, 
                                  outputFastq = paste0(result_dir_parsed,  "/fastp_dual/", file_pair$common_part, "/processed_", file_pair$common_part),
                                  thread = 2,
                                  verbose = TRUE
  )
  return(rfastp_result)
}



# Get values and dataframes
extractValues <- function(results){
  summary = results %>% lapply(extract_summary_parts)
  filtering_result = results %>% 
    lapply(function(x){x$filtering_result}) %>% 
    lapply(as.data.frame) %>% 
    lapply(t) %>% 
    lapply(as.data.frame) %>% 
    lapply(rownames_to_column) 
  duplication = results %>% lapply(function(x){x$duplication})
  insert_size = results %>% lapply(function(x){x$insert_size})
  adapter_cutting = results %>% 
    lapply(function(x){x$adapter_cutting}) %>%
    lapply(as.data.frame) %>% 
    lapply(t) %>% 
    lapply(as.data.frame) %>% 
    lapply(rownames_to_column) 
  read1_before_filtering = results %>% lapply(function(x){x$read1_before_filtering}) %>% lapply(extract_whatever_filtering)
  read2_before_filtering = results %>% lapply(function(x){x$read2_before_filtering}) %>% lapply(extract_whatever_filtering)
  read1_after_filtering = results %>% lapply(function(x){x$read1_after_filtering}) %>% lapply(extract_whatever_filtering)
  read2_after_filtering = results %>% lapply(function(x){x$read2_after_filtering})%>% lapply(extract_whatever_filtering)
  
  result_list = list(summary, filtering_result, duplication, insert_size, adapter_cutting, read1_before_filtering, read2_before_filtering, read1_after_filtering, read2_after_filtering)
  return(result_list)
}


