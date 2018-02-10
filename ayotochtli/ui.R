#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

source("helper.r")


# Define UI for application
load("data.Rdata")  # for the dataset

shinyUI(pageWithSidebar(
  headerPanel('Data'),
  sidebarPanel(
     
    radioButtons('x', "X-axis",
                 list("ANOVA FDR" = 9,
                      "ANOVA pval" = 10,
                      "Mouse DEG recurrence" = 11,
                      "Human DEG recurrence" = 12)),

    radioButtons('logx', "log?",
                 list("default" = 0, 
                      "-log10" = 1,
                      "log10" = 2)),
    
    radioButtons('y', "Y-axis",
                 list("ANOVA FDR" = 9,
                      "ANOVA pval" = 10,
                      "Mouse DEG recurrence" = 11,
                      "Human DEG recurrence" = 12)),
    
    radioButtons('logy', "log?",
                 list("default" = 0, 
                      "-log10" = 1,
                      "log10" = 2)),
    
    sliderInput("bins",
                "Number of bins:",
                min = 1,
                max = 50,
                value = 30),
    br(),
    checkboxGroupInput('show_vars', 'Columns to show:', names(gene_table),
                       selected = names(gene_table))
 
  ),
  
  
  mainPanel(
    tabsetPanel(
      tabPanel("Plot", plotOutput("myplot")), 
      tabPanel("BeanPlots", plotOutput("mybarplot")), 
      tabPanel("Histograms", plotOutput("myhist")), 
      tabPanel('Genes of interest', dataTableOutput("mytable"))
    )
  )
))



