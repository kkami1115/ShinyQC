#install.packages("pacman")
#pacman::p_load(BiocManager, dplyr)
#BiocManager::install(c("Rfastp"))
#library(Rfastp)

# functions

#' Create file pairs list
#'
#' @param paths A character vector of file paths
#' @param pattern A regular expression pattern to match file names
#' @return A list of file pairs
#' @importFrom magrittr %>%
create_file_pairs_list <- function(paths, pattern) {
  common_parts <- basename(gsub(paste0(pattern, "\\.fastq\\.gz"), "", paths))
  unique_common_parts <- unique(common_parts)

  file_pairs_list <- list()
  for (common_part in unique_common_parts) {
    matched_files <- grep(common_part, paths, value = TRUE)
    if (length(matched_files) == 2) {
      file_pairs_list[[common_part]] <- list(
        common_part = common_part,
        file1 = matched_files[1],
        file2 = matched_files[2]
      )
    }
  }

  return(file_pairs_list)
}

#' Extract summary parts
#'
#' @param youso A list containing summary information
#' @return A data frame with 'before_filtering' and 'after_filtering' columns
extract_summary_parts <- function(youso) {
  before_filtering_df <- youso$summary$before_filtering %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column()

  after_filtering_df <- youso$summary$after_filtering %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column()

  colnames(before_filtering_df) <- c("item", "before_filtering")
  colnames(after_filtering_df) <- c("item", "after_filtering")

  dplyr::inner_join(before_filtering_df, after_filtering_df, by = "item")
}

#' Extract filtering information
#'
#' @param youso A list containing filtering information
#' @return A list with filtering details
extract_filtering_info <- function(youso) {
  list(
    main = data.frame(
      total_reads = youso$total_reads,
      total_bases = youso$total_bases,
      q20_bases = youso$q20_bases,
      q30_bases = youso$q30_bases
    ),
    quality_curves = as.data.frame(youso$quality_curves),
    content_curves = as.data.frame(youso$content_curves),
    kmer_count = t(as.data.frame(youso$kmer_count))
  )
}


#' Run Rfastp
#'
#' @param result_dir_parsed Parsed result directory path
#' @param file_pair A list containing file pair information
#' @param threads Number of threads to use (default: 2)
#' @return Rfastp result object
run_rfastp <- function(result_dir_parsed, file_pair, threads = 2) {
  output_dir <- file.path(result_dir_parsed, "fastp_dual", file_pair$common_part)
  dir.create(output_dir, recursive = TRUE, mode = "0777")

  output_fastq <- file.path(output_dir, paste0("trimmed_", file_pair$common_part))

  Rfastp::rfastp(
    read1 = file_pair$file1,
    read2 = file_pair$file2,
    outputFastq = output_fastq,
    thread = threads,
    verbose = TRUE
  )
}



#' Extract values and data frames from Rfastp results
#'
#' @param results A list of Rfastp result objects
#' @return A list containing extracted values and data frames
extract_values <- function(results) {
  list(
    summary = lapply(results, extract_summary_parts),
    filtering_result = lapply(results, function(x) {
      x$filtering_result %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column()
    }),
    duplication = lapply(results, function(x) x$duplication),
    insert_size = lapply(results, function(x) x$insert_size),
    adapter_cutting = lapply(results, function(x) {
      x$adapter_cutting %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column()
    }),
    read1_before_filtering = lapply(results, function(x) extract_filtering_info(x$read1_before_filtering)),
    read2_before_filtering = lapply(results, function(x) extract_filtering_info(x$read2_before_filtering)),
    read1_after_filtering = lapply(results, function(x) extract_filtering_info(x$read1_after_filtering)),
    read2_after_filtering = lapply(results, function(x) extract_filtering_info(x$read2_after_filtering))
  )
}
