library(shiny)
source('fitting.R')

# Define UI for application that draws a histogram
ui <- navbarPage("Flood Frequency Analysis",
     tabPanel("Explore",
     
       sidebarLayout(
          sidebarPanel(
            # choose gauging station
            selectInput("station", "Station",
                        c("Jackson" = "USGS02486000", 
                          "Edinburg" = "USGS02482000", 
                          "Carthage" = "USGS02482550", 
                          "Lena" = "USGS02483500", 
                          "Rockport" = "USGS02488000", 
                          "Monticello" = "USGS02488500", 
                          "Columbia" = "USGS02489000",
                          "Bogalusa" = "USGS02489500")),
            
            # select plot type
            selectInput("plottype", "Select",
                        c("Time series", 
                          "Density"))
          ),
          
          mainPanel(
             plotOutput("timePlot")
          )
        )
     ),
     tabPanel("Fit",
          sidebarLayout(
            sidebarPanel(
              # choose distribution
              selectInput("distr", "Distribution",
                          c("Normal" = "norm", 
                            "Log-normal" = "lognorm",
                            "Exponential" = "expo",
                            "Gamma" = "gam",
                            "Pearson 3" = "p3",
                            "Log-Pearson 3" = "lp3",
                            "Gumbel" = "gum",
                            "Weibull" = "wei")),
              
              # choose estimation method
              radioButtons("method",
                      "Estimation method",
                      c("Maximum Likelihood (MLE)" = "mle", 
                        "Method of Moments (MOM)" = "mme", 
                        "Probablity weighted moments (PWM)" = "pwme")),
              h5("RMSE:"),
              uiOutput('rmse'),
              h5("K-S:"),
              uiOutput('ks'),
              h5("A-D:"),
              uiOutput('ad'),
              h5("AIC:"),
              uiOutput('aic'),
              h5("BIC:"),
              uiOutput('bic'),
              h5("Parameters"),
              verbatimTextOutput('params')
              
            ),
            mainPanel(
              plotOutput("fitPlot"),
              plotOutput("quantilePlot"),
              plotOutput("finalPlot") 
            )
          )
     )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    dataInput <- reactive({ 
      data <- read_excel(paste0("data/", input$station, ".xls"), sheet = 1) 
      out <- preprocessing(data)
      data_processed <<- out$processed
      data_sorted <<- out$sorted
      aps <<- out$aps
      return_periods <<- (dim(data_sorted)[1]+1)/data_sorted$rank    # Weibull plotting position   
      if (input$distr == "norm") {
        estim <<- normal(aps, input$method, log = FALSE, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "lognorm") {
        estim <<- normal(log(aps), input$method, log = TRUE, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "expo") {
        estim <<- exponential(aps, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "gam") {
        estim <<- gam(aps, input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "p3") {
        estim <<- p3(aps, input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "lp3") {
        estim <<- lp3(log(aps), input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "gum") {
        estim <<- gumb(aps, input$method, rp = return_periods, plotx = xseq)
      }
      if (input$distr == "wei") {
        estim <<- weibull(aps, input$method, rp = return_periods, plotx = xseq)
      }
      })
    
    output$timePlot <- renderPlot({
      dataInput()
      if (input$plottype == "Time series") {
          # time series plot
          plot(data_processed$Year, data_processed$aps, type="l", xlab = "Year", ylab = "Maximum discharge (cf/s)", main = "Time series of the raw data")
      }
      
      if (input$plottype == "Density") {  
          hist(aps, breaks=15, main = "Distribution of annual peak discharge", xlab = "Annual peak streamflow (cf/s)", prob = TRUE)
          lines(density(aps), col = "blue", lwd = 2)
      }
  })

    output$fitPlot <- renderPlot({
      dataInput()
      hist(aps, breaks=15, main = "Distribution of annual peak discharge", xlab = "Annual peak streamflow (cf/s)", prob = TRUE)
      lines(density(aps), col = "blue", lwd = 2)
      lines(xseq, estim$ploty, col = "red" , type="l", lwd = 2)
    })
    
    output$quantilePlot <- renderPlot({
      dataInput()
      qqplot(aps, estim$xt, xlab="Observed peaks (cfs)", ylab="Estimated peaks (cfs)", main="Q-Q plot")
      abline(c(0,1), col="red")
    })
    
    output$finalPlot <- renderPlot({
      dataInput()
      plot(log(return_periods), data_sorted$aps, xlab = "Log return period", ylab = "Discharge (cfs)", 
           main = "Estimated peak discharge vs log return period")
      lines(log(return_periods), estim$xt, col=2)
    })
    
    output$ks <- renderPrint({
      dataInput()
      cat(round(ks.test(estim$xt, aps)$statistic, digits=3))
    })
    
    output$rmse <- renderPrint({
      dataInput()
      cat(round(sqrt(sum((sort(aps, decreasing = TRUE)-estim$xt)^2)/length(aps)), digits=1))
    })
    
    output$aic <- renderPrint({
      dataInput()
      if(input$method != "mle"){
        cat("Only available for MLE") 
      } else {
      cat(round(-2*estim$likelihood+2*estim$parno, digits=1))
      }
    })
    
    output$bic <- renderPrint({
      dataInput()
      if(input$method != "mle"){
       cat("Only available for MLE")
      } else {
        cat(round(-2*estim$likelihood+estim$parno*log(length(aps)), digits=1))
      }
    })
    
    output$ad <- renderPrint({
      dataInput()
      cat(round(estim$ad, digits=2))
    })
    
    output$params <- renderPrint({
      dataInput()
      names <- names(estim$par)
      for (i in 1:estim$parno){
        cat(names[i], ": ", toString(signif(unlist(estim$par)[i], 5)))
        cat("\n")
      }
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

