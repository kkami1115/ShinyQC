
pacman::p_load(shiny, shinyFiles, shinyjs)
if (require(bslib)) {
  theme <- bslib::bs_theme(version = 4)
} else {
  theme <- NULL
}

useShinyjs()

fluidPage(
  
  titlePanel("Fastp execution and visualization"),
  tabsetPanel(id="tabs",
    tabPanel("Fastp execution",
             # fastqフォルダ選択用ボタン
             shinyDirButton("fastq_dir", "Fastq Directory", "Select fastq directory"),
             verbatimTextOutput("fastq_dir"),
             shinyDirButton("result_dir", "Result Directory", "Select result directory"),
             verbatimTextOutput("result_dir"),
             textInput("file_pattern", "File Pattern", value = "_[12]" ),
             actionButton("runButton", "Run Rfastp as paired")
    ),
    tabPanel("Overview",
             selectInput("overview_item", "Select item", choices = c("total_reads", "total_bases", "q20_bases", "q30_bases", "q20_rate", "q30_rate", "read1_mean_length", "read2_mean_length", "gc_content") ),
             plotlyOutput("all_samples_summary")
    ),
    tabPanel("Visualize per sample", 
             sidebarLayout(
               sidebarPanel(
                 selectInput("selected_sample", "Select sample:", choices = NULL)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary", tableOutput("summary")),
                   tabPanel("Filtering Results", tableOutput("filtering_result_table")),
                   tabPanel("Duplicate Rate",
                            verbatimTextOutput("duplication_rate"),
                            plotlyOutput("duplication_histogram"),
                            plotlyOutput("duplication_mean_gc")
                   ),
                   tabPanel("Insert Length",
                            verbatimTextOutput("insert_size_peak"),
                            verbatimTextOutput("insert_size_unknown"),
                            plotlyOutput("insert_size_histogram")
                   ),
                   tabPanel("Adapter cutting",
                            tableOutput("adapter_cutting")        
                   ),
                   tabPanel("Read1 Before Filtering",
                            # read1 before filtering関連の表示
                            tableOutput("read1_before_filtering_main"),
                            plotlyOutput("read1_before_filtering_qualitycurves"),
                            plotlyOutput("read1_before_filtering_contentcurves"),
                            plotlyOutput("read1_before_filtering_kmer")
                   ),
                   tabPanel("Read1 After Filtering",
                            # read1 after filtering関連の表示
                            tableOutput("read1_after_filtering_main"),
                            plotlyOutput("read1_after_filtering_qualitycurves"),
                            plotlyOutput("read1_after_filtering_contentcurves"),
                            plotlyOutput("read1_after_filtering_kmer")
                   ),
                   tabPanel("Read2 Before Filtering",
                            # read2 before filtering関連の表示
                            tableOutput("read2_before_filtering_main"),
                            plotlyOutput("read2_before_filtering_qualitycurves"),
                            plotlyOutput("read2_before_filtering_contentcurves"),
                            plotlyOutput("read2_before_filtering_kmer")
                   ),
                   tabPanel("Read2 After Filtering",
                            # read2 after filtering関連の表示
                            tableOutput("read2_after_filtering_main"),
                            plotlyOutput("read2_after_filtering_qualitycurves"),
                            plotlyOutput("read2_after_filtering_contentcurves"),
                            plotlyOutput("read2_after_filtering_kmer")
                   )
                 )
               )
             )
    )
  )
)

