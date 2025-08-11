
#pacman::p_load(shiny, shinyFiles, shinyjs, plotly, bslib)
#source("./R/fastp.R")

theme <- bslib::bs_theme(version = 4)

shiny::fluidPage(
  shinyjs::useShinyjs(),
  shiny::titlePanel("Fastp execution and visualization"),

  # tabset panel
  shiny::tabsetPanel(id="tabs",

    # Fastp execution tab
    shiny::tabPanel("Fastp execution",
                   # Button for selecting fastq dir
                   shinyFiles::shinyDirButton("fastq_dir", "Fastq Directory", "Select a directory containing fastq files"),
                   shiny::verbatimTextOutput("fastq_dir"),
                   shinyFiles::shinyDirButton("result_dir", "Result Directory", "Select a directory to save the results"),
                   shiny::verbatimTextOutput("result_dir"),
                   shiny::textInput("file_pattern", "File Pattern", value = "_[12]" ),
                   shiny::actionButton("runButton", "Run Rfastp (Paired-end Mode)")
                   #,
                   # Area of logs
                   # uiOutput("dynamicSizeLogOutput")
    ),

    # Overview tab
    shiny::tabPanel("Overview",
                    shiny::sidebarLayout(
                      shiny::sidebarPanel(
                        shiny::selectInput("overview_item", "Select Quality Metric",
                                           choices = c("total_reads",
                                                       "total_bases",
                                                       "q20_bases",
                                                       "q30_bases",
                                                       "q20_rate",
                                                       "q30_rate",
                                                       "read1_mean_length",
                                                       "read2_mean_length",
                                                       "gc_content")
                        )
                      ),
                      shiny::mainPanel(
                        plotly::plotlyOutput("all_samples_summary")
                      )
                    )
    ),


    # visualize per sample tab


    shiny::tabPanel("Visualize per sample",
                    shiny::sidebarLayout(
                      shiny::sidebarPanel(
                        shiny::selectInput("selected_sample", "Select sample:", choices = NULL)
                      ),
                      shiny::mainPanel(
                        shiny::tabsetPanel(
                          shiny::tabPanel("Summary", shiny::tableOutput("summary")),
                          shiny::tabPanel("Filtering Results", shiny::tableOutput("filtering_result_table")),
                          shiny::tabPanel("Duplicate Rate",
                                          shiny::verbatimTextOutput("duplication_rate"),
                                          plotly::plotlyOutput("duplication_histogram"),
                                          plotly::plotlyOutput("duplication_mean_gc")
                          ),
                          shiny::tabPanel("Insert Length",
                                          shiny::verbatimTextOutput("insert_size_peak"),
                                          shiny::verbatimTextOutput("insert_size_unknown"),
                                          plotly::plotlyOutput("insert_size_histogram")
                          ),
                          shiny::tabPanel("Adapter cutting",
                                          shiny::tableOutput("adapter_cutting")
                          ),
                          shiny::tabPanel("Read1 Before Filtering",
                                          shiny::tableOutput("read1_before_filtering_main"),
                                          plotly::plotlyOutput("read1_before_filtering_qualitycurves"),
                                          plotly::plotlyOutput("read1_before_filtering_contentcurves"),
                                          plotly::plotlyOutput("read1_before_filtering_kmer")
                          ),
                          shiny::tabPanel("Read1 After Filtering",
                                          shiny::tableOutput("read1_after_filtering_main"),
                                          plotly::plotlyOutput("read1_after_filtering_qualitycurves"),
                                          plotly::plotlyOutput("read1_after_filtering_contentcurves"),
                                          plotly::plotlyOutput("read1_after_filtering_kmer")
                          ),
                          shiny::tabPanel("Read2 Before Filtering",
                                          shiny::tableOutput("read2_before_filtering_main"),
                                          plotly::plotlyOutput("read2_before_filtering_qualitycurves"),
                                          plotly::plotlyOutput("read2_before_filtering_contentcurves"),
                                          plotly::plotlyOutput("read2_before_filtering_kmer")
                          ),
                          shiny::tabPanel("Read2 After Filtering",
                                          shiny::tableOutput("read2_after_filtering_main"),
                                          plotly::plotlyOutput("read2_after_filtering_qualitycurves"),
                                          plotly::plotlyOutput("read2_after_filtering_contentcurves"),
                                          plotly::plotlyOutput("read2_after_filtering_kmer")
                          )
                        )
                      )
                    )
    ),

    # Settings tab
    # Settings tab
    shiny::tabPanel("Settings",
                    shiny::selectInput("multi_strategy", "Choose a strategy:",
                                       list("sequential", "multisession")),
                    shiny::selectInput("threads", "Thread(s) for fastp",
                                       list(1,2,3,4,5,6,7,8))
                      )
  )
)

