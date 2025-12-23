# install.packages("pacman")
# pacman::p_load(BiocManager, dplyr)
# BiocManager::install(c("Rfastp"))
# library(Rfastp)

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

# create_file_pairs_list <- function(paths, pattern) {
#   df <- tibble(path = paths,
#                common_part = basename(gsub(paste0(pattern, "\\.fastq\\.gz"), "", paths)))
#   pairs <- df %>% group_by(common_part) %>% filter(n() == 2) %>%
#     summarise(file1 = first(path), file2 = last(path), .groups = "drop")
#   setNames(
#     purrr::pmap(list(pairs$common_part, pairs$file1, pairs$file2),
#                 ~list(common_part = ..1, file1 = ..2, file2 = ..3)),
#     pairs$common_part
#   )
# }


#' Fix integer overflow in summary dataframe
#'
#' @param df A dataframe with columns 'item' and 'value'
#' @param is_before Boolean, true if this is before_filtering data
#' @param target_mean Numeric, target mean length (used for after_filtering)
#' @return A dataframe with corrected values
fix_overflow_df <- function(df, is_before = TRUE, target_mean = NULL) {
  # Helper to get value
  get_val <- function(k) as.numeric(df[df$item == k, 2])
  set_val <- function(k, v) df[df$item == k, 2] <<- v

  total_reads <- get_val("total_reads")
  total_bases <- get_val("total_bases")

  # Fix total_bases
  # Initial correction to positive range [0, 2^32] if negative
  if (total_bases < 0) total_bases <- total_bases + 2^32

  current_mean <- total_bases / total_reads

  if (is_before) {
    # Match against standard NGS read lengths
    standards <- c(36, 50, 75, 100, 125, 150, 200, 250, 300)
    best_k <- 0
    min_diff <- Inf

    # Try adding k * 2^32 (k=0, 1, 2...)
    for (k in 0:5) {
      test_bases <- total_bases + k * 2^32
      test_mean <- test_bases / total_reads

      # Find distance to nearest standard
      diff <- min(abs(test_mean - standards))

      if (diff < min_diff) {
        min_diff <- diff
        best_k <- k
      }
    }
    total_bases <- total_bases + best_k * 2^32
  } else {
    # For after_filtering, aim for mean length <= target_mean but close to it
    # Trimming shouldn't drastically reduce length usually, but can.
    # We assume mean_length_after <= mean_length_before

    if (!is.null(target_mean)) {
      best_k <- 0
      min_diff <- Inf

      for (k in 0:5) {
        test_bases <- total_bases + k * 2^32
        test_mean <- test_bases / total_reads

        if (test_mean <= target_mean + 5) { # Allow slight tolerance
          diff <- target_mean - test_mean
          if (diff < min_diff) {
            min_diff <- diff
            best_k <- k
          }
        }
      }
      total_bases <- total_bases + best_k * 2^32
    }
  }
  set_val("total_bases", total_bases)

  # Fix q20, q30 bases
  # They should be <= total_bases.
  # If negative, add 2^32.
  # Then add 2^32 until reasonable ratio?
  # Assumption: Q20/Q30 rate is usually high (>50%).
  # If val / total_bases < 0.3 (and we added 2^32), maybe add another?
  # But simpler: ensure val <= total_bases.
  # And maximize val given val <= total_bases?

  for (k in c("q20_bases", "q30_bases")) {
    val <- get_val(k)
    if (val < 0) val <- val + 2^32

    # Try adding 2^32 as long as it stays <= total_bases
    # This maximizes the Q-score bases count within the total_bases limit
    while ((val + 2^32) <= total_bases) {
      val <- val + 2^32
    }
    set_val(k, val)
  }

  # Recalculate rates
  q20_bases <- get_val("q20_bases")
  q30_bases <- get_val("q30_bases")
  set_val("q20_rate", q20_bases / total_bases)
  set_val("q30_rate", q30_bases / total_bases)

  # Fix mean lengths if negative or suspicious
  r1_len <- get_val("read1_mean_length")
  r2_len <- get_val("read2_mean_length")

  # If we fixed total_bases, we should probably update mean lengths to be consistent
  # especially if they were negative.
  # Even if positive, they might be wrong if they overflowed differently?
  # But usually we trust our total_bases fix more.
  avg_len <- total_bases / total_reads

  if (r1_len < 0 || r2_len < 0) {
    set_val("read1_mean_length", avg_len)
    set_val("read2_mean_length", avg_len)
  } else {
    # If reported length is very different from calculated, maybe update?
    # Let's update if difference is large (> 10bp)
    if (abs(r1_len - avg_len) > 10) set_val("read1_mean_length", avg_len)
    if (abs(r2_len - avg_len) > 10) set_val("read2_mean_length", avg_len)
  }

  return(df)
}

#' Extract summary parts
#'
#' @param youso A list containing summary information
#' @param fix_overflow Boolean, whether to apply overflow fix (default: TRUE)
#' @return A data frame with 'before_filtering' and 'after_filtering' columns
extract_summary_parts <- function(youso, fix_overflow = TRUE) {
  before_filtering_df <- youso$summary$before_filtering %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column()

  simple_fix_overflow <- function(df) {
    target_items <- c("total_reads", "total_bases", "q20_bases", "q30_bases")
    for (item in target_items) {
      val <- as.numeric(df[df$item == item, 2])
      if (!is.na(val) && val < 0) {
        df[df$item == item, 2] <- val + 2^32
      }
    }

    # Recalculate rates
    total_bases <- as.numeric(df[df$item == "total_bases", 2])
    q20_bases <- as.numeric(df[df$item == "q20_bases", 2])
    q30_bases <- as.numeric(df[df$item == "q30_bases", 2])

    if (!is.na(total_bases) && total_bases > 0) {
      df[df$item == "q20_rate", 2] <- q20_bases / total_bases
      df[df$item == "q30_rate", 2] <- q30_bases / total_bases
    }

    return(df)
  }

  colnames(before_filtering_df) <- c("item", "value")
  if (fix_overflow) {
    before_filtering_df <- fix_overflow_df(before_filtering_df, is_before = TRUE)
  } else {
    before_filtering_df <- simple_fix_overflow(before_filtering_df)
  }

  # Get target mean from fixed before df
  tr_b <- as.numeric(before_filtering_df[before_filtering_df$item == "total_reads", 2])
  tb_b <- as.numeric(before_filtering_df[before_filtering_df$item == "total_bases", 2])
  target_mean <- tb_b / tr_b

  colnames(before_filtering_df) <- c("item", "before_filtering")

  after_filtering_df <- youso$summary$after_filtering %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column()

  colnames(after_filtering_df) <- c("item", "value")
  if (fix_overflow) {
    after_filtering_df <- fix_overflow_df(after_filtering_df, is_before = FALSE, target_mean = target_mean)
  } else {
    after_filtering_df <- simple_fix_overflow(after_filtering_df)
  }
  colnames(after_filtering_df) <- c("item", "after_filtering")

  dplyr::inner_join(before_filtering_df, after_filtering_df, by = "item")
}


#' Extract filtering information
#'
#' @param youso A list containing filtering information
#' @return A list with filtering details
extract_filtering_info <- function(youso) {
  fix_neg <- function(x) {
    if (is.numeric(x) && x < 0) x + 2^32 else x
  }

  process_curves <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      as.data.frame(lapply(x, unlist))
    } else {
      as.data.frame(x)
    }
  }

  process_kmer <- function(x) {
    if (is.list(x) && !is.data.frame(x)) {
      df <- as.data.frame(unlist(x))
      colnames(df) <- "V1"
      df
    } else {
      t(as.data.frame(x))
    }
  }

  list(
    main = data.frame(
      total_reads = fix_neg(youso$total_reads),
      total_bases = fix_neg(youso$total_bases),
      q20_bases = fix_neg(youso$q20_bases),
      q30_bases = fix_neg(youso$q30_bases)
    ),
    quality_curves = process_curves(youso$quality_curves),
    content_curves = process_curves(youso$content_curves),
    kmer_count = process_kmer(youso$kmer_count)
  )
}


#' Run Rfastp
#'
#' @param output_dir  output directory
#' @param file_pair A list containing file pair information
#' @param threads Number of threads to use (default: 2)
#' @return Rfastp result object
run_rfastp <- function(output_dir, file_pair, threads = 2) {
  dir.create(file.path(output_dir, file_pair$common_part), recursive = TRUE, mode = "0777")

  output_fastq <- file.path(output_dir, file_pair$common_part, paste0("trimmed_", file_pair$common_part))

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
#' @param fix_overflow Boolean, whether to apply overflow fix (default: TRUE)
#' @return A list containing extracted values and data frames
extract_values <- function(results, fix_overflow = TRUE) {
  list(
    summary = lapply(results, function(x) extract_summary_parts(x, fix_overflow = fix_overflow)),
    filtering_result = lapply(results, function(x) {
      x$filtering_result %>%
        as.data.frame() %>%
        t() %>%
        as.data.frame() %>%
        tibble::rownames_to_column() %>%
        mutate(V1 = as.numeric(V1))
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
# extract_values <- function(results) {
#   list(
#     summary = purrr::map(results, extract_summary_parts),
#     filtering_result = purrr::map(results,
#                                   ~ .x$filtering_result %>% as.data.frame() %>% t() %>% as.data.frame() %>%
#                                     tibble::rownames_to_column() %>% dplyr::mutate(V1 = as.numeric(V1))),
#     duplication = purrr::map(results, "duplication"),
#     insert_size = purrr::map(results, "insert_size"),
#     adapter_cutting = purrr::map(results,
#                                  ~ .x$adapter_cutting %>% as.data.frame() %>% t() %>% as.data.frame() %>%
#                                    tibble::rownames_to_column()),
#     read1_before_filtering = purrr::map(results, ~ extract_filtering_info(.x$read1_before_filtering)),
#     read2_before_filtering = purrr::map(results, ~ extract_filtering_info(.x$read2_before_filtering)),
#     read1_after_filtering  = purrr::map(results, ~ extract_filtering_info(.x$read1_after_filtering)),
#     read2_after_filtering  = purrr::map(results, ~ extract_filtering_info(.x$read2_after_filtering))
#   )
# }



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

# convert_units_summary <- function(df) {
#   df %>%
#     mutate(
#       before_filtering = case_when(
#         item %in% c("total_reads","total_bases","q20_bases","q30_bases") ~ convert_value_with_unit(before_filtering),
#         item %in% c("q20_rate","q30_rate","gc_content") ~ paste0(round(before_filtering * 100, 2), "%"),
#         item %in% c("read1_mean_length","read2_mean_length") ~ paste0(round(before_filtering), " bp"),
#         TRUE ~ as.character(before_filtering)
#       ),
#       after_filtering = case_when(
#         item %in% c("total_reads","total_bases","q20_bases","q30_bases") ~ convert_value_with_unit(after_filtering),
#         item %in% c("q20_rate","q30_rate","gc_content") ~ paste0(round(after_filtering * 100, 2), "%"),
#         item %in% c("read1_mean_length","read2_mean_length") ~ paste0(round(after_filtering), " bp"),
#         TRUE ~ as.character(after_filtering)
#       )
#     )
# }


#' convert digits for filtering_result
#'
#' @param df filtering_result df
#' @return digits-modified filtering_result df
convert_units_filtering_result <- function(df) {
  passed_filter_reads <- round(df[df$rowname == "passed_filter_reads", "V1"])
  low_quality_reads <- round(df[df$rowname == "low_quality_reads", "V1"])
  too_many_N_reads <- round(df[df$rowname == "too_many_N_reads", "V1"])
  too_short_reads <- round(df[df$rowname == "too_short_reads", "V1"])
  too_long_reads <- round(df[df$rowname == "too_long_reads", "V1"])

  df[df$item == "passed_filter_reads", 2] <- convert_value_with_unit(passed_filter_reads)
  df[df$item == "low_quality_reads", 2] <- convert_value_with_unit(low_quality_reads)
  df[df$item == "too_many_N_reads", 2] <- convert_value_with_unit(too_many_N_reads)
  df[df$item == "too_short_reads", 2] <- convert_value_with_unit(too_short_reads)
  df[df$item == "too_long_reads", 2] <- round(too_long_reads, 0)

  colnames(df) <- c("item", "value")
  return(df)
}

# convert_units_filtering_result <- function(df) {
#   df %>%
#     mutate(
#       value = case_when(
#         rowname %in% c("passed_filter_reads","low_quality_reads","too_many_N_reads","too_short_reads") ~ convert_value_with_unit(round(V1)),
#         rowname == "too_long_reads" ~ as.character(round(V1,0)),
#         TRUE ~ as.character(V1)
#       )
#     ) %>%
#     rename(item = rowname) %>%
#     select(item, value)
# }
