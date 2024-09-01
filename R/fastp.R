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
#' @param output_dir  output directory
#' @param file_pair A list containing file pair information
#' @param threads Number of threads to use (default: 2)
#' @return Rfastp result object
run_rfastp <- function(result_dir_parsed, output_dir, file_pair, threads = 2) {
  dir.create( file.path(output_dir, file_pair$common_part), recursive = TRUE, mode = "0777")

  output_fastq <- file.path(output_dir, file_pair$common_part,  paste0("trimmed_", file_pair$common_part))

  result <- Rfastp::rfastp(
    read1 = file_pair$file1,
    read2 = file_pair$file2,
    outputFastq = output_fastq,
    thread = threads,
    verbose = TRUE
  )

  return(result)
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

convert_value_with_unit <- function(value) {
  if (is.numeric(value)) {
    if (abs(value) >= 1e9) {
      return(paste0(round(value / 1e9, 2), " G"))
    } else if (abs(value) >= 1e6) {
      return(paste0(round(value / 1e6, 2), " M"))
    } else {
      return(as.character(value))
    }
  } else {
    return(as.character(value))
  }
}

#' convert digits for summary
#'
#' @param df summary df
#' @return digits-modified summary df
convert_units_summary <- function(df) {
  total_reads_before <- df[df$item == "total_reads", "before_filtering"]
  total_reads_after <- df[df$item == "total_reads", "after_filtering"]

  total_bases_before <- df[df$item == "total_bases", "before_filtering"]
  total_bases_after <- df[df$item == "total_bases", "after_filtering"]

  q20_bases_before <- df[df$item == "q20_bases", "before_filtering"]
  q20_bases_after <- df[df$item == "q20_bases", "after_filtering"]

  q30_bases_before <- df[df$item == "q30_bases", "before_filtering"]
  q30_bases_after <- df[df$item == "q30_bases", "after_filtering"]

  q20_rate_before <- df[df$item == "q20_rate", "before_filtering"] * 100
  q20_rate_after <- df[df$item == "q20_rate", "after_filtering"] * 100

  q30_rate_before <- df[df$item == "q30_rate", "before_filtering"] * 100
  q30_rate_after <- df[df$item == "q30_rate", "after_filtering"] * 100

  read1_length_before <- round(df[df$item == "read1_mean_length", "before_filtering"])
  read1_length_after <- round(df[df$item == "read1_mean_length", "after_filtering"])

  read2_length_before <- round(df[df$item == "read2_mean_length", "before_filtering"])
  read2_length_after <- round(df[df$item == "read2_mean_length", "after_filtering"])

  gc_content_before <- df[df$item == "gc_content", "before_filtering"] * 100
  gc_content_after <- df[df$item == "gc_content", "after_filtering"] * 100

  df[df$item == "total_reads", 2:3] <- c(convert_value_with_unit(total_reads_before), convert_value_with_unit(total_reads_after))
  df[df$item == "total_bases", 2:3] <- c(convert_value_with_unit(total_bases_before), convert_value_with_unit(total_bases_after))
  df[df$item == "q20_bases", 2:3] <- c(convert_value_with_unit(q20_bases_before), convert_value_with_unit(q20_bases_after))
  df[df$item == "q30_bases", 2:3] <- c(convert_value_with_unit(q30_bases_before), convert_value_with_unit(q30_bases_after))
  df[df$item == "q20_rate", 2:3] <- c(paste0(round(q20_rate_before, 2), "%"), paste0(round(q20_rate_after, 2), "%"))
  df[df$item == "q30_rate", 2:3] <- c(paste0(round(q30_rate_before, 2), "%"), paste0(round(q30_rate_after, 2), "%"))
  df[df$item == "read1_mean_length", 2:3] <- c(paste0(read1_length_before, " bp"), paste0(read1_length_after, " bp"))
  df[df$item == "read2_mean_length", 2:3] <- c(paste0(read2_length_before, " bp"), paste0(read2_length_after, " bp"))
  df[df$item == "gc_content", 2:3] <- c(paste0(round(gc_content_before, 2), "%"), paste0(round(gc_content_after, 2), "%"))

  return(df)
}


#' convert digits for filtering_result
#'
#' @param df filtering_result df
#' @return digits-modified filtering_result df
convert_units_filtering_result <- function(df) {
  passed_filter_reads <- df[df$rowname == "passed_filter_reads", "V1"]
  low_quality_reads <- df[df$rowname == "low_quality_reads", "V1"]
  too_many_N_reads <- df[df$rowname == "too_many_N_reads", "V1"]
  too_short_reads	 <- df[df$rowname == "too_short_reads", "V1"]
  too_long_reads <- df[df$rowname == "too_long_reads", "V1"]

  df[df$item == "passed_filter_reads", 2] <- convert_value_with_unit(passed_filter_reads)
  df[df$item == "low_quality_reads", 2] <- convert_value_with_unit(low_quality_reads)
  df[df$item == "too_many_N_reads", 2] <- convert_value_with_unit(too_many_N_reads)
  df[df$item == "too_short_reads", 2] <- convert_value_with_unit(too_short_reads)
  df[df$item == "too_long_reads", 2] <- round(too_long_reads, 0)

  colnames(df)  = c("item", "value")
  return(df)
}

