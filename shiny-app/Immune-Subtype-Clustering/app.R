#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
source('src/load_deps.R')

options(shiny.maxRequestSize=50*1024^2)

imagePlot <- function(m) {
  print(m)
  print(colnames(m))
  image(1:ncol(m), 1:nrow(m), t(m), col = terrain.colors(60), 
        axes = FALSE, main='PanCan Cluster Label Check (percent samples)', 
        xlab = 'Reported Cluster Labels', ylab='New Cluster Calls')
  axis(1, 1:ncol(m), colnames(m))
  axis(2, 1:nrow(m), rownames(m))
  for (x in 1:ncol(m))
    for (y in 1:nrow(m))
      text(x, y, m[y,x])
}

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("PanImmune Subtypes"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      # Input: Select separator ----
      radioButtons("sep", "File Separator",
                   choices = c(Comma = ",",
                               Tab = "\t"),
                   selected = ","),
      checkboxInput("logged", "Apply Log10", TRUE),
      
      # Horizontal line ----
      tags$hr(),
      
      # Input: Select a file ----
      fileInput("file1", "Choose CSV file with \n1st column gene symbols",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",
                           ".csv.gz",
                           "text/tsv",
                           "text/comma-separated-values,text/plain",
                           ".tsv",
                           ".tsv.gz"),
                placeholder = 'data/ivy20.csv'),
      
      # Horizontal line ----
      tags$hr(),
      
      numericInput("corenum", "Cores", 4, width = '100'),
      
      actionButton("gobutton", "GO"),
      h4("WARNING! This is in beta! \n Median centered RSEM RPKM gene expression values were used to train the model.")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      # table output for the scores and subtype for each input sample
      # maybe... pick a sample, and get the score x subtype distr plot?
      # distribution of input data compared to TCGA data
      plotOutput("distPlot"),
      plotOutput('barPlot'),
      fluidRow(
        DT::dataTableOutput("table")
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # in src files ... have same path as app.R
  reportedClusters <- getTable()
  
  # get new calls
  getCalls <- eventReactive(input$gobutton, {
    newdat <- input$file1
    withProgress(message = 'Working...', value = 0, {
      newScores(newdat, input$logged, input$corenum)
    })
  })
  
  # plot of where a sample is in signature space X clusters    
  output$distPlot <- renderPlot({
    #heatmap(as.matrix(getCalls()$Table), xlab = 'Reported Clusters', ylab = 'New Calls')
    imagePlot(getCalls()$Table)
  })
  
  output$barPlot <- renderPlot({
    counts <- table(getCalls()$MaxCalls)
    barplot(counts, main="New Cluster Label Calls", 
            xlab="Cluster Labels")
  })
  
  
  # Filter data based on selections
  output$table <- DT::renderDataTable(
    DT::datatable(
      as.data.frame(getCalls()$ProbCalls),
      extensions = 'Buttons', options = list(
        dom = 'Bfrtip',
        buttons = 
          list('copy', 'print', 
               list(
                 extend = 'collection',
                 buttons = c('csv', 'excel', 'pdf'),
                 text = 'Download')
          )
      )
      
    )
  )
}



bioc <- local({
  env <- new.env()
  on.exit(rm(env))
  evalq(source("http://bioconductor.org/biocLite.R", local = TRUE), env)
  biocinstallRepos()
})

r <- getOption("repos")
for (i in 1:length(bioc)) {
  print(bioc)
  r[names(bioc)[i]] <- bioc[i]
}
options(repos = as.character(r))
options(repos = BiocInstaller::biocinstallRepos())

source('src/tablef_fun.R')
source('src/computing_scores_and_calling_clusters.R')


# Run the application 
shinyApp(ui = ui, server = server)



