
pacman::p_load(shiny,shinyjs,fs, future,furrr)
source("fastp.R")

# 非同期処理の計画を設定
plan(multisession)  # 複数のプロセスを使用


shinyServer(function(input, output, session) {
  
  ################## ファイル指示系 ############################
  # フォルダ選択
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

  # 表示系
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
  
  
  ################## Rfastp実行系 ############################  
  # 
  observeEvent(input$runButton, {
    shinyjs::disable("tabs")
    fastq_dir_parsed <- parseDirPath(volumes,input$fastq_dir)
    result_dir_parsed <- parseDirPath(volumes, input$result_dir)
    
    # リザルトディレクトリを作成
    dir.create( paste0(result_dir, "/fastp_dual/"), showWarnings = FALSE, recursive = TRUE  )
    
    # Fastqファイルのパスを取得
    fastq_paths <- list.files(path = fastq_dir_parsed, pattern = "\\.fastq\\.gz$", full.names = TRUE, recursive = TRUE)
    # ファイルペアを生成
    file_pairs_list <- createFilePairsList(fastq_paths, pattern = input$file_pattern)
    
    # 結果を格納するリスト
    results <- list()
    
    # Rfastpの実行
    results <- future_map(file_pairs_list, ~runRfastp(.x)) 
    result_list <- extractValues(results)
    shinyjs::enable()
  })
  
  
  
  ################# Overview以降可視化部分 #################
  summary= result_list[[1]]
  filtering_result = result_list[[2]]
  duplication = result_list[[3]] 
  insert_size = result_list[[4]]
  adapter_cutting = result_list[[5]]
  read1_before_filtering = result_list[[6]]
  read2_before_filtering = result_list[[7]]
  read1_after_filtering = result_list[[8]]
  read2_after_filtering = result_list[[9]]
  
  
  
  output$all_samples_summary <- renderPlot({
    tmp = lapply(names(summary), function(sample_name) {
      df <- summary[[sample_name]]
      df <- df[df$item == input$overview_item, ]
      df$sample <- sample_name
      return(df)
    }) %>% dplyr::bind_rows() %>%
      pivot_longer(cols = c("before_filtering", "after_filtering"))
    tmp$name =  factor(tmp$name, levels = c("before_filtering", "after_filtering"))
    
    ggplot(tmp, aes(x=sample, y=value, fill=name)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_fill_manual(values = c("before_filtering" = "blue", "after_filtering" = "red")) 
  })
  
  output$summary <- renderTable({
    summary[[input$selected_sample]]
  })
  output$filtering_result_table <- renderTable({
    filtering_result[[input$selected_sample]]
  })
  
  
  # duplicate関連の表示
  output$duplication_rate <- renderText({
    rate <- duplication[[input$selected_sample]]$rate
    paste("duplicate率:", rate)
  })
  output$duplication_histogram <- renderPlot({
    histogram_data <- duplication[[input$selected_sample]]$histogram
    barplot(histogram_data, main="Histogram of Duplication", xlab="Duplication", ylab="Frequency")
  })
  output$duplication_mean_gc <- renderPlot({
    mean_gc_data <- duplication[[input$selected_sample]]$mean_gc
    plot(mean_gc_data, type='o', col='blue', main="Mean GC Content", xlab="Position", ylab="Mean GC")
  })
  
  
  # insert size関連の表示
  output$insert_size_peak <- renderText({
    peak <- insert_size[[input$selected_sample]]$peak
    paste("insert size peak:", peak)
  })
  output$insert_size_unknown <- renderText({
    unknown <- insert_size[[input$selected_sample]]$unknown
    paste("insert size unknown:", unknown)
  })
  output$insert_size_histogram <- renderPlot({
    insert_size_histogram <- insert_size[[input$selected_sample]]$histogram
    plot(insert_size_histogram, type='o', col='blue', main="Insert size", xlab="Position", ylab="value")
  })
  
  
  # adapter cutting関連の表示
  output$adapter_cutting <- renderTable({
    adapter_cutting[[input$selected_sample]] 
  })
  
  
  # read1 before filtering関連の表示
  output$read1_before_filtering_main <- renderTable({
    read1_before_filtering[[input$selected_sample]]$main
  })
  output$read1_before_filtering_qualitycurves <- renderPlot({
    tmp <- read1_before_filtering[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read1_before_filtering_contentcurves <- renderPlot({
    tmp <- read1_before_filtering[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read1_before_filtering_kmer <- renderPlot({
    tmp <- read1_before_filtering[[input$selected_sample]]$kmer_count
    tmp %>%  as.data.frame() %>%
      rownames_to_column() %>% 
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>% 
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>% 
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile() 
  })
  
  # read1 after filtering関連の表示
  output$read1_after_filtering_main <- renderTable({
    read1_after_filtering[[input$selected_sample]]$main
  })
  output$read1_after_filtering_qualitycurves <- renderPlot({
    tmp <- read1_after_filtering[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read1_after_filtering_contentcurves <- renderPlot({
    tmp <- read1_after_filtering[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read1_after_filtering_kmer <- renderPlot({
    tmp <- read1_after_filtering[[input$selected_sample]]$kmer_count
    tmp %>%  as.data.frame() %>%
      rownames_to_column() %>% 
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>% 
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>% 
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile() 
  })
  
  # read2 before filtering関連の表示
  output$read2_before_filtering_main <- renderTable({
    read2_before_filtering[[input$selected_sample]]$main
  })
  output$read2_before_filtering_qualitycurves <- renderPlot({
    tmp <- read2_before_filtering[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read2_before_filtering_contentcurves <- renderPlot({
    tmp <- read2_before_filtering[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read2_before_filtering_kmer <- renderPlot({
    tmp <- read2_before_filtering[[input$selected_sample]]$kmer_count
    tmp %>%  as.data.frame() %>%
      rownames_to_column() %>% 
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>% 
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>% 
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile() 
  })
  
  # read2 after filtering関連の表示
  output$read2_after_filtering_main <- renderTable({
    read2_after_filtering[[input$selected_sample]]$main
  })
  output$read2_after_filtering_qualitycurves <- renderPlot({
    tmp <- read2_after_filtering[[input$selected_sample]]$quality_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read2_after_filtering_contentcurves <- renderPlot({
    tmp <- read2_after_filtering[[input$selected_sample]]$content_curves  %>% pivot_longer(cols = everything())
    ggplot(tmp, aes(x=1:length(name), y=value, color=name)) + geom_line() 
  })
  output$read2_after_filtering_kmer <- renderPlot({
    tmp <- read2_after_filtering[[input$selected_sample]]$kmer_count
    tmp %>%  as.data.frame() %>%
      rownames_to_column() %>% 
      mutate(first_char = str_sub(.data$rowname, 1 ,2)) %>% 
      mutate(sore_igai = str_sub(.data$rowname, start = 3) ) %>% 
      ggplot(aes(x = first_char, y = sore_igai, fill = V1)) +
      geom_tile() 
  })
})
