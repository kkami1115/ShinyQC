
#pacman::p_load(shiny,shinyjs,fs, future,furrr, plotly, progressr)
source("fastp.R")
source("utils-pipe.R")
source("library_tmp.R")

# Specify seed is fixed
.options = furrr_options(seed = TRUE)

# Set handlers for progress bar
progressr::handlers(progressr::handler_shiny)


shiny::shinyServer(function(input, output, session) {


  ################## Directories section ############################
  # Choose the directories
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())
  shinyFiles::shinyDirChoose(input, "fastq_dir", roots = volumes, session = session)
  shinyFiles::shinyDirChoose(input, "result_dir", roots = volumes, session = session)

  # display
  output$fastq_dir <- shiny::renderPrint({
    if (is.integer(input$fastq_dir)) {
      cat("No directory have been selected (fastq_dir)")
    } else {
      shinyFiles::parseDirPath(volumes, input$fastq_dir)
    }
  })

  output$result_dir <- shiny::renderPrint({
    if (is.integer(input$result_dir)) {
      cat("No directory has been selected (result_dir)")
    } else {
      shinyFiles::parseDirPath(volumes, input$result_dir)
    }
  })


  ################## Rfastp exec section ############################
  # handles results value dynamically
  result_values <- shiny::reactiveValues(result_list = NULL)


  shiny::observeEvent(input$runButton, {
    shinyjs::disable("tabs")

    threads = as.numeric(input$threads)

    # Set plan for asynchronous processing
    if(input$multi_strategy == "sequential"){
      future::plan(future::sequential())   # Process one file at a time
    }else if(input$multi_strategy == "multisession"){
      future::plan(future::multisession(), workers = (parallel::detectCores() / threads))  # MaxThreads/2 since multiple sessions are launched and rfastp runs on 2 threads.
    }

    # Set result directories
    fastq_dir_parsed <- shinyFiles::parseDirPath(volumes,input$fastq_dir)
    result_dir_parsed <- shinyFiles::parseDirPath(volumes, input$result_dir)
    output_dir <- file.path(result_dir_parsed, paste("fastp_dual",format(Sys.Date(), "%Y%m%d"), sep = "_") )

    # Validate input
    if (is.null(fastq_dir_parsed) || is.null(result_dir_parsed)) {
      shiny::showNotification("Please select a valid directory.", type = "error", closeButton = TRUE)
      return()
    }
    # Get fastq files fullpath and check existence
    fastq_paths <- list.files(fastq_dir_parsed, pattern = "\\.fastq\\.gz$|\\.fastq$|\\.fq$", full.names = TRUE, recursive = TRUE)
    if (length(fastq_paths) == 0) {
      shiny::showNotification("We can't find any fastq files in the directory.", type = "error", closeButton = TRUE)
    }


    # Make fastq filepairs as R{1,2}
    file_pairs_list <- create_file_pairs_list(fastq_paths, pattern = input$file_pattern)
    if (length(file_pairs_list) == 0) {
      shiny::showNotification("File pairing failed or no valid file pairs found.", type = "error", closeButton = TRUE)
    }

    # Dynamically update `selectInput` choices
    shiny::observe({
      shiny::req(file_pairs_list) # Check `result_summary` isn't NULL
      shiny::updateSelectInput(session, "selected_sample",
                               choices = names(file_pairs_list))
    })

    # Popup while processing
    shiny::showModal(shiny::modalDialog(
      title = "Processing",
      "Please wait...",
      footer = NULL
    ))

    # list for results
    results <- list()

    # function for exec run_rfastp
    rfastp_exec <- function(p, x, threads){
      furrr::future_map(x, ~{
        p(.x$common_part)
        run_rfastp(result_dir_parsed, output_dir, .x, threads)
      },
      seed = TRUE, globals = TRUE)
    }


    # exec above function
    progressr::withProgressShiny(message = "Processing...",
       {
         p <- progressr::progressor(along = names(file_pairs_list))
         results <- rfastp_exec(p, file_pairs_list, threads)
       }
    )

    result_list <- extract_values(results)
    result_values$result_list <- result_list

    # enable tabs
    shinyjs::enable("tabs")

    if (!is.null(result_list)){
      # Notification for complete
      shiny::showModal(shiny::modalDialog(
        title = "Complete",
        "Rfastp processing is complete.",
        easyClose = TRUE,
        footer = shiny::modalButton("Close")
      ))
    }

  })



  ################# Visualization section #################
  result_summary= shiny::reactive({result_values$result_list[[1]]})
  filtering_result = shiny::reactive({result_values$result_list[[2]]})
  duplication = shiny::reactive({result_values$result_list[[3]] })
  insert_size = shiny::reactive({result_values$result_list[[4]]})
  adapter_cutting = shiny::reactive({result_values$result_list[[5]]})
  read1_before_filtering = shiny::reactive({result_values$result_list[[6]]})
  read2_before_filtering = shiny::reactive({result_values$result_list[[7]]})
  read1_after_filtering = shiny::reactive({result_values$result_list[[8]]})
  read2_after_filtering = shiny::reactive({result_values$result_list[[9]]})


  output$all_samples_summary <- plotly::renderPlotly({
    summary_data <- result_summary()
    req(summary_data)
    tmp <- lapply(names(summary_data), function(sample_name) {
      df <- summary_data[[sample_name]]
      df <- df[df$item == input$overview_item, ]
      df$sample <- sample_name
      return(df)
    }) %>% dplyr::bind_rows() %>%
      tidyr::pivot_longer(cols = c("before_filtering", "after_filtering"))
    tmp$name =  factor(tmp$name, levels = c("before_filtering", "after_filtering"))

    p <- ggplot2::ggplot(tmp, aes(x=sample, y=value, fill=name)) +
      ggplot2::geom_bar(stat = "identity", position = position_dodge()) +
      ggplot2::scale_fill_manual(values = c("before_filtering" = "blue", "after_filtering" = "red"))
    p %>% plotly::ggplotly()
  })


  output$summary <- shiny::renderTable({
    summary_data <- result_summary()
    req(summary_data)
    summary_data[[input$selected_sample]] %>% convert_units_summary()
  })

  output$filtering_result_table <- shiny::renderTable({
    filtering_data <- filtering_result()
    req(filtering_data)
    filtering_data[[input$selected_sample]]  %>% lapply(., function(col) {sapply(col, convert_value_with_unit) }) %>% as.data.frame()
  })


  # elements of duplicate rates
  output$duplication_rate <- shiny::renderText({
    rate <- duplication()[[input$selected_sample]]$rate
    req(rate)
    paste("duplicate rate:", round(rate * 100 , digits = 2), " %")
  })

  output$duplication_histogram <- plotly::renderPlotly({
    histogram_data <- duplication()[[input$selected_sample]]$histogram
    # barplot(histogram_data, xlab="Duplication Level", ylab="Read percent (%)")
    req(histogram_data)
    p <- histogram_data %>%
      as.data.frame() %>%
      tibble::rowid_to_column() %>%
      dplyr::rename(c("x"="rowid","y"=".")) %>%
      ggplot2::ggplot(aes(x=x,y=y, group=x)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::xlab("Duplication Level") +
      ggplot2::ylab("Read Count (raw)")
    p %>% plotly::ggplotly()
  })

  output$duplication_mean_gc <- plotly::renderPlotly({
    mean_gc_data <- duplication()[[input$selected_sample]]$mean_gc
    # plot(mean_gc_data, type='o', col='blue', main="Mean GC Content", xlab="Position", ylab="Mean GC")
    req(mean_gc_data)
    p <- mean_gc_data %>%
      as.data.frame() %>%
      mutate(across(everything(), ~ . * 100 )) %>%
      tibble::rowid_to_column() %>%
      dplyr::rename(c("x"="rowid","y"=".")) %>%
      ggplot2::ggplot(aes(x=x,y=y)) +
      ggplot2::geom_line() +
      ggplot2::xlab("Duplication Level") +
      ggplot2::ylab("Mean GC Ratio (%)")
    p %>% plotly::ggplotly()
  })


  # elements of insert size
  output$insert_size_peak <- shiny::renderText({
    peak <- insert_size()[[input$selected_sample]]$peak
    req(peak)
    paste("insert size peak:", peak)
  })
  output$insert_size_unknown <- shiny::renderText({
    unknown <- insert_size()[[input$selected_sample]]$unknown
    req(unknown)
    paste("insert size unknown:", unknown %>% convert_value_with_unit() )
  })
  output$insert_size_histogram <- plotly::renderPlotly({
    insert_size_histogram <- insert_size()[[input$selected_sample]]$histogram
    # plot(insert_size_histogram, type='o', col='blue', main="Insert size", xlab="Position", ylab="value")
    req(insert_size_histogram)
    p <- insert_size_histogram %>%
      as.data.frame() %>%
      tibble::rowid_to_column() %>%
      dplyr::rename(c("x"="rowid","y"=".")) %>%
      ggplot2::ggplot(aes(x=x,y=y, group=x)) + ggplot2::geom_bar(stat = "identity")
    p %>% plotly::ggplotly()
  })


  # elements of adapter cutting
  output$adapter_cutting <- shiny::renderTable({
    req(adapter_cutting()[[input$selected_sample]])
    adapter_cutting()[[input$selected_sample]] %>% lapply(., function(col) {sapply(col, convert_value_with_unit) }) %>% as.data.frame()
  })


  # elements of read1 before filtering
  output$read1_before_filtering_main <- shiny::renderTable({
    read1_before_filtering()[[input$selected_sample]]$main %>% lapply(., function(col) {sapply(col, convert_value_with_unit) }) %>% as.data.frame()
  })
  output$read1_before_filtering_qualitycurves <- plotly::renderPlotly({
    req(read1_before_filtering()[[input$selected_sample]]$quality_curves )
        tmp <- read1_before_filtering()[[input$selected_sample]]$quality_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_before_filtering_contentcurves <- plotly::renderPlotly({
    req(read1_before_filtering()[[input$selected_sample]]$content_curves )
    tmp <- read1_before_filtering()[[input$selected_sample]]$content_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_before_filtering_kmer <- plotly::renderPlotly({
    req(read1_before_filtering()[[input$selected_sample]]$kmer_count )
    tmp <- read1_before_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(char_1st2nd = stringr::str_sub(.data$rowname, 1 ,2)) %>%
      dplyr::mutate(others = stringr::str_sub(.data$rowname, start = 3) ) %>%
      ggplot2::ggplot(aes(x = char_1st2nd, y = others, fill = V1)) +
      ggplot2::geom_tile()
    p %>% plotly::ggplotly()
  })

  # elements of read1 after filtering
  output$read1_after_filtering_main <- shiny::renderTable({
    read1_after_filtering()[[input$selected_sample]]$main %>% lapply(., function(col) {sapply(col, convert_value_with_unit) }) %>% as.data.frame()
  })
  output$read1_after_filtering_qualitycurves <- plotly::renderPlotly({
    req(read1_after_filtering()[[input$selected_sample]]$quality_curves)
    tmp <- read1_after_filtering()[[input$selected_sample]]$quality_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_after_filtering_contentcurves <- plotly::renderPlotly({
    req(read1_after_filtering()[[input$selected_sample]]$content_curves)
    tmp <- read1_after_filtering()[[input$selected_sample]]$content_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_after_filtering_kmer <- plotly::renderPlotly({
    req(read1_after_filtering()[[input$selected_sample]]$kmer_count)
    tmp <- read1_after_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(char_1st2nd = stringr::str_sub(.data$rowname, 1 ,2)) %>%
      dplyr::mutate(others = stringr::str_sub(.data$rowname, start = 3) ) %>%
      ggplot2::ggplot(aes(x = char_1st2nd, y = others, fill = V1)) +
      ggplot2::geom_tile()
    p %>% plotly::ggplotly()
  })

  # elements of read2 before filtering
  output$read2_before_filtering_main <- shiny::renderTable({
    read2_before_filtering()[[input$selected_sample]]$main %>% lapply(., function(col) {sapply(col, convert_value_with_unit) }) %>% as.data.frame()
  })
  output$read2_before_filtering_qualitycurves <- plotly::renderPlotly({
    req(read2_before_filtering()[[input$selected_sample]]$quality_curves)
    tmp <- read2_before_filtering()[[input$selected_sample]]$quality_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_before_filtering_contentcurves <- plotly::renderPlotly({
    req(read2_before_filtering()[[input$selected_sample]]$content_curves)
    tmp <- read2_before_filtering()[[input$selected_sample]]$content_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_before_filtering_kmer <- plotly::renderPlotly({
    req(read2_before_filtering()[[input$selected_sample]]$kmer_count)
    tmp <- read2_before_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(char_1st2nd = stringr::str_sub(.data$rowname, 1 ,2)) %>%
      dplyr::mutate(others = stringr::str_sub(.data$rowname, start = 3) ) %>%
      ggplot2::ggplot(aes(x = char_1st2nd, y = others, fill = V1)) +
      ggplot2::geom_tile()
    p %>% plotly::ggplotly()
  })

  # elements of read2 after filtering
  output$read2_after_filtering_main <- shiny::renderTable({
    read2_after_filtering()[[input$selected_sample]]$main %>% lapply(., function(col) {sapply(col, convert_value_with_unit) }) %>% as.data.frame()
  })
  output$read2_after_filtering_qualitycurves <- plotly::renderPlotly({
    req(read2_after_filtering()[[input$selected_sample]]$quality_curves)
    tmp <- read2_after_filtering()[[input$selected_sample]]$quality_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) +
      ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_after_filtering_contentcurves <- plotly::renderPlotly({
    req(read2_after_filtering()[[input$selected_sample]]$content_curves)
    tmp <- read2_after_filtering()[[input$selected_sample]]$content_curves  %>%
      tidyr::pivot_longer(cols = everything())
    p <- ggplot2::ggplot(tmp, aes(x=1:length(name), y=value, color=name)) +
      ggplot2::geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_after_filtering_kmer <- plotly::renderPlotly({
    req(read2_after_filtering()[[input$selected_sample]]$kmer_count)
    tmp <- read2_after_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      dplyr::mutate(char_1st2nd = stringr::str_sub(.data$rowname, 1 ,2)) %>%
      dplyr::mutate(others = stringr::str_sub(.data$rowname, start = 3) ) %>%
      ggplot2::ggplot(aes(x = char_1st2nd, y = others, fill = V1)) +
      ggplot2::geom_tile()
    p %>% plotly::ggplotly()
  })

  # Back future::plan to sequential (default)
  future::plan(future::sequential())
})
