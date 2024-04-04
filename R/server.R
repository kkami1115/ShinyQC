
#pacman::p_load(shiny,shinyjs,fs, future,furrr, plotly, progressr)
source("./R/fastp.R")

# Set handlers for progress bar
progressr::handlers(progressr::handler_shiny)


shiny::shinyServer(function(input, output, session) {


  ################## Folders set section ############################
  # Choose the directories
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())
  shinyFiles::shinyDirChoose(input, "fastq_dir", roots = volumes, session = session)
  shinyFiles::shinyDirChoose(input, "result_dir", roots = volumes, session = session)

  observe({
    cat("\ninput$fastq_dir value:\n\n")
    print(input$fastq_dir)
  })
  observe({
    cat("\ninput$result_dir value:\n\n")
    print(input$result_dir)
  })

  # display
  output$fastq_dir <- renderPrint({
    if (is.integer(input$fastq_dir)) {
      cat("No files have been selected (fastq_dir)")
    } else {
      parseDirPath(volumes, input$fastq_dir)
    }
  })

  output$result_dir <- renderPrint({
    if (is.integer(input$result_dir)) {
      cat("No directory has been selected (result_dir)")
    } else {
      parseDirPath(volumes, input$result_dir)
    }
  })


  ################## Rfastp exec section ############################
  # handles results value dynamically
  result_values <- reactiveValues(result_list = NULL)
  # handles logfiles path dynamically
  # log_paths <- reactiveValues(files= NULL)

  observeEvent(input$runButton, {
    shinyjs::disable("tabs")
    #showNotification("Processing has started...", type = "message", duration = 5) # display for 5sec

    # Set plan for asynchronous processing
    if(input$multi_strategy == "sequential"){
      plan(sequential)   # Process one file at a time
    }else if(input$multi_strategy == "multisession"){
      plan(multisession, workers = parallel::detectCores() / 2)  # MaxThreads/2 since multiple sessions are launched and rfastp runs on 2 threads.
    }

    # Set result directories
    fastq_dir_parsed <- parseDirPath(volumes,input$fastq_dir)
    result_dir_parsed <- parseDirPath(volumes, input$result_dir)

    # Make directories
    dir.create( paste0(result_dir_parsed, "/fastp_dual/"), showWarnings = FALSE, recursive = TRUE  )

    # Get fastq files fullpath
    fastq_paths <- list.files(path = fastq_dir_parsed, pattern = "\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE)
    # Make fastq filepairs as R{1,2}
    file_pairs_list <- createFilePairsList(fastq_paths, pattern = input$file_pattern)

    # Dynamically update `selectInput` choices
    observe({
      req(file_pairs_list) # Check `result_summary` isn't NULL
      updateSelectInput(session, "selected_sample",
                        choices = names(file_pairs_list))
    })

    # Popup while processing
    showModal(modalDialog(
     title = "Processing",
     "Please wait...",
     footer = NULL
    ))

    # list for results
    results <- list()

    # function for exec runRfastp
    rfastp_exec <- function(p, x){
     furrr::future_map(x, ~{
        p(.x$common_part)
        runRfastp(result_dir_parsed, .x)
        },
        seed = 1L, globals = TRUE)
    }

    # exec above function
    withProgressShiny(message = "Now calculating...",
      {
       p <- progressor(along = names(file_pairs_list))
       results <- rfastp_exec(p, file_pairs_list)
    }
    )

    # results <- future_map(file_pairs_list, ~runRfastp(result_dir_parsed, .x), seed=1L, globals=TRUE )
    result_list <- extractValues(results)
    result_values$result_list <- result_list

    # re-enable tabs
    shinyjs::enable("tabs")

    # Notification for complete
    showModal(modalDialog(
      title = "Complete",
      "Rfastp processing is complete.",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))

  })



  ################# Visualization section #################
  result_summary= reactive({result_values$result_list[[1]]})
  filtering_result = reactive({result_values$result_list[[2]]})
  duplication = reactive({result_values$result_list[[3]] })
  insert_size = reactive({result_values$result_list[[4]]})
  adapter_cutting = reactive({result_values$result_list[[5]]})
  read1_before_filtering = reactive({result_values$result_list[[6]]})
  read2_before_filtering = reactive({result_values$result_list[[7]]})
  read1_after_filtering = reactive({result_values$result_list[[8]]})
  read2_after_filtering = reactive({result_values$result_list[[9]]})



  output$all_samples_summary <- renderPlotly({
    summary_data <- result_summary()
    tmp <- lapply(names(summary_data), function(sample_name) {
      df <- summary_data[[sample_name]]
      df <- df[df$item == input$overview_item, ]
      df$sample <- sample_name
      return(df)
    }) %>% dplyr::bind_rows() %>%
      pivot_longer(cols = c("before_filtering", "after_filtering"))
    tmp$name =  factor(tmp$name, levels = c("before_filtering", "after_filtering"))

    p <- ggplot(tmp, aes(x=sample, y=value, fill=name)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_fill_manual(values = c("before_filtering" = "blue", "after_filtering" = "red"))
    p %>% plotly::ggplotly()
  })



  output$summary <- renderTable({
    summary_data <- result_summary()
    # req(!is.null(summary_data), nrow(summary_data) > 0)
    summary_data[[input$selected_sample]]
  })

  output$filtering_result_table <- renderTable({
    filtering_data <- filtering_result()
    # req(!is.null(filtering_data), nrow(filtering_data) > 0)
    filtering_data[[input$selected_sample]]
  })


  # elements of duplicate rates
  output$duplication_rate <- renderText({
    rate <- duplication()[[input$selected_sample]]$rate
    paste("duplicate rate:", rate)
  })
  output$duplication_histogram <- renderPlotly({
    histogram_data <- duplication()[[input$selected_sample]]$histogram
    # barplot(histogram_data, xlab="Duplication Level", ylab="Read percent (%)")
    p <- histogram_data %>%
      as.data.frame() %>%
      rowid_to_column() %>%
      rename(c("x"="rowid","y"=".")) %>%
      ggplot(aes(x=x,y=y, group=x)) +
      geom_bar(stat = "identity") +
      xlab("Duplication Level") +
      ylab("Read Count (raw)")
    p %>% ggplotly()
  })
  output$duplication_mean_gc <- renderPlotly({
    mean_gc_data <- duplication()[[input$selected_sample]]$mean_gc
    # plot(mean_gc_data, type='o', col='blue', main="Mean GC Content", xlab="Position", ylab="Mean GC")
    p <- mean_gc_data %>%
      as.data.frame() %>%
      rowid_to_column() %>%
      rename(c("x"="rowid","y"=".")) %>%
      ggplot(aes(x=x,y=y)) +
      geom_line() +
      xlab("Duplication Level") +
      ylab("Mean GC Ratio (%)")
    p %>% ggplotly()
  })


  # elements of insert size
  output$insert_size_peak <- renderText({
    peak <- insert_size()[[input$selected_sample]]$peak
    paste("insert size peak:", peak)
  })
  output$insert_size_unknown <- renderText({
    unknown <- insert_size()[[input$selected_sample]]$unknown
    paste("insert size unknown:", unknown)
  })
  output$insert_size_histogram <- renderPlotly({
    insert_size_histogram <- insert_size()[[input$selected_sample]]$histogram
    # plot(insert_size_histogram, type='o', col='blue', main="Insert size", xlab="Position", ylab="value")
    p <- insert_size_histogram %>%
      as.data.frame() %>%
      rowid_to_column() %>%
      rename(c("x"="rowid","y"=".")) %>%
      ggplot(aes(x=x,y=y, group=x)) + geom_bar(stat = "identity")
    p %>% plotly::ggplotly()
  })


  # elements of adapter cutting
  output$adapter_cutting <- renderTable({
    adapter_cutting()[[input$selected_sample]]
  })


  # elements of read1 before filtering
  output$read1_before_filtering_main <- renderTable({
    read1_before_filtering()[[input$selected_sample]]$main
  })
  output$read1_before_filtering_qualitycurves <- renderPlotly({
    tmp <- read1_before_filtering()[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_before_filtering_contentcurves <- renderPlotly({
    tmp <- read1_before_filtering()[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_before_filtering_kmer <- renderPlotly({
    tmp <- read1_before_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>%
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>%
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile()
    p %>% plotly::ggplotly()
  })

  # elements of read1 after filtering
  output$read1_after_filtering_main <- renderTable({
    read1_after_filtering()[[input$selected_sample]]$main
  })
  output$read1_after_filtering_qualitycurves <- renderPlotly({
    tmp <- read1_after_filtering()[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_after_filtering_contentcurves <- renderPlotly({
    tmp <- read1_after_filtering()[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read1_after_filtering_kmer <- renderPlotly({
    tmp <- read1_after_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>%
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>%
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile()
    p %>% plotly::ggplotly()
  })

  # elements of read2 before filtering
  output$read2_before_filtering_main <- renderTable({
    read2_before_filtering()[[input$selected_sample]]$main
  })
  output$read2_before_filtering_qualitycurves <- renderPlotly({
    tmp <- read2_before_filtering()[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_before_filtering_contentcurves <- renderPlotly({
    tmp <- read2_before_filtering()[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_before_filtering_kmer <- renderPlotly({
    tmp <- read2_before_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>%
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>%
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile()
    p %>% plotly::ggplotly()
  })

  # elements of read2 after filtering
  output$read2_after_filtering_main <- renderTable({
    read2_after_filtering()[[input$selected_sample]]$main
  })
  output$read2_after_filtering_qualitycurves <- renderPlotly({
    tmp <- read2_after_filtering()[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_after_filtering_contentcurves <- renderPlotly({
    tmp <- read2_after_filtering()[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    p <- ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line()
    p %>% plotly::ggplotly()
  })
  output$read2_after_filtering_kmer <- renderPlotly({
    tmp <- read2_after_filtering()[[input$selected_sample]]$kmer_count
    p <- tmp %>%  as.data.frame() %>%
      tibble::rownames_to_column() %>%
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>%
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>%
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile()
    p %>% plotly::ggplotly()
  })

  # Back plan to sequential (default)
  plan(sequential)
})
