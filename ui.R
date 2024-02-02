
pacman::p_load(shiny, shinyFiles, shinyjs)
if (require(bslib)) {
  theme <- bslib::bs_theme(version = 4)
} else {
  theme <- NULL
}


fluidPage(
  
  titlePanel("Fastp execution and visualization"),
  tabsetPanel(id="tabs",
    tabPanel("Fastp execution",
             # fastqフォルダ選択用ボタン
             shinyDirButton("fastq_dir", "Fastq Directory", "Fastqファイルが存在するフォルダを選択"),
             verbatimTextOutput("fastq_dir"),
             shinyDirButton("result_dir", "Result Directory", "結果を保存するフォルダを選択"),
             verbatimTextOutput("result_dir"),
             textInput("file_pattern", "File Pattern", value = "_[12]" ),
             actionButton("runButton", "Run Rfastp as paired")
    ),
    tabPanel("Overview",
             selectInput("overview_item", "アイテムを選択", choices = c("total_reads", "total_bases", "q20_bases", "q30_bases", "q20_rate", "q30_rate", "read1_mean_length", "read2_mean_length", "gc_content") ),
             plotOutput("all_samples_summary")
    ),
    tabPanel("Visualize per sample", 
             sidebarLayout(
               sidebarPanel(
                 selectInput("selected_sample", "サンプルを選択:", choices = names(summary))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput("summary")),
                   tabPanel("Filtering Results", tableOutput("filtering_result_table")),
                   tabPanel("Duplicate Rate",
                            verbatimTextOutput("duplication_rate"),
                            plotOutput("duplication_histogram"),
                            plotOutput("duplication_mean_gc")
                   ),
                   tabPanel("Insert Length",
                            verbatimTextOutput("insert_size_peak"),
                            verbatimTextOutput("insert_size_unknown"),
                            plotOutput("insert_size_histogram")
                   ),
                   tabPanel("Adapter cutting",
                            tableOutput("adapter_cutting")        
                   ),
                   tabPanel("Read1 Before Filtering",
                            # read1 before filtering関連の表示
                            tableOutput("read1_before_filtering_main"),
                            plotOutput("read1_before_filtering_qualitycurves"),
                            plotOutput("read1_before_filtering_contentcurves"),
                            plotOutput("read1_before_filtering_kmer")
                   ),
                   tabPanel("Read1 After Filtering",
                            # read1 after filtering関連の表示
                            tableOutput("read1_after_filtering_main"),
                            plotOutput("read1_after_filtering_qualitycurves"),
                            plotOutput("read1_after_filtering_contentcurves"),
                            plotOutput("read1_after_filtering_kmer")
                   ),
                   tabPanel("Read2 Before Filtering",
                            # read2 before filtering関連の表示
                            tableOutput("read2_before_filtering_main"),
                            plotOutput("read2_before_filtering_qualitycurves"),
                            plotOutput("read2_before_filtering_contentcurves"),
                            plotOutput("read2_before_filtering_kmer")
                   ),
                   tabPanel("Read2 After Filtering",
                            # read2 after filtering関連の表示
                            tableOutput("read2_after_filtering_main"),
                            plotOutput("read2_after_filtering_qualitycurves"),
                            plotOutput("read2_after_filtering_contentcurves"),
                            plotOutput("read2_after_filtering_kmer")
                   )
                 )
               )
             )
    )
  )
)

