#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

shinyServer(function(input, output) {
  load("data.Rdata")  
  
  datax <- reactive({  
    switch(input$x,
                   "9" = 9,
                   "10" = 10,
                   "11" = 11,
                   "12" = 12)
  })    
 datay <- reactive({  
      
    switch(input$y,
                "9" = 9,
                "10" = 10,
                "11" = 11,
                "12" = 12)
    
  })
 
 logdatax <- reactive({  
   switch(input$logx,
          "0" = 0,
          "1" = 1,
          "2" = 2)
 }) 
 
 logdatay <- reactive({  
   
   switch(input$logy,
          "0" = 0,
          "1" = 1,
          "2" = 2)
   
 })
 
 
 
  output$myplot <- renderPlot({
    x = datax()
    y = datay()
    logx = logdatax()
    logy = logdatay()
    dataX = gene_table[,x]
    dataY = gene_table[,y] 
    if(logx==1){  dataX = -log10(dataX) }
    if(logx==2){  dataX = log10(dataX) }
    if(logy==1){  dataY = -log10(dataY) }
    if(logy==2){  dataY = log10(dataY) }
    a=names(gene_table)[x]
    b=names(gene_table)[y]
    
    plot(dataX, dataY, pch=19, bty="n", xlab=a, ylab=b)

  })
  
  output$myhist <- renderPlot({
    x = datax()
    logx = logdatax()
    dataX = gene_table[,x]
    if(logx==1){  dataX = -log10(dataX) }
    if(logx==2){  dataX = log10(dataX) }
    a=names(gene_table)[x]
    dataX[!is.finite(dataX)] = NA
    # generate bins based on input$bins from ui.R
    bins <- seq(min(dataX, na.rm=T), max(dataX, na.rm=T), length.out = input$bins + 1)
    
    # draw the histogram with the specified number of bins
    hist(dataX, breaks = bins, col = viridis(input$bins+1), border = 'white', main="", xlab=a)    
  })
  
  output$mybarplot <- renderPlot({
    x = datax()
    y = datay()
    logx = logdatax()
    logy = logdatay()
    dataX = gene_table[,x]
    dataY = gene_table[,y] 
    if(logx==1){  dataX = -log10(dataX) }
    if(logx==2){  dataX = log10(dataX) }
    if(logy==1){  dataY = -log10(dataY) }
    if(logy==2){  dataY = log10(dataY) }
    
    a=names(gene_table)[x]
    b=names(gene_table)[y]
    
    beanplot(list(dataX, dataY), 
             what=c(1,1,1,0), col=list("green", "purple"), bty="n", names=list(a,b))
    
  })
  # a large table, reative to input$show_vars
  output$mytable = renderDataTable({

    gene_table[, input$show_vars, drop = FALSE]
  })
  
  
})