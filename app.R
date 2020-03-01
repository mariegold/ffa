library(shiny)
source('fitting.R')

# Define UI for application that draws a histogram
ui <- navbarPage("Flood Frequency Analysis",
     tabPanel("Explore",
     
       sidebarLayout(
          sidebarPanel(
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
                        "Probablity weighted moments (PWM)" = "pwme"))
            ),
            mainPanel(
              plotOutput("fitPlot"),
              plotOutput("quantilePlot") 
              #           hover = "plot_hover"),
              # verbatimTextOutput("info")
            )
          )
     )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$timePlot <- renderPlot({
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
      hist(aps, breaks=15, main = "Distribution of annual peak discharge", xlab = "Annual peak streamflow (cf/s)", prob = TRUE)
      
      if (input$distr == "norm") {  
        est <- normal(aps, input$method, log = FALSE)
        lines(xseq, dnorm(xseq, est$mu, est$sigma), col = "red", type="l", lwd = 2)
      }
      if (input$distr == "lognorm") {  
        est <- normal(aps, input$method, log = TRUE)
        lines(xseq, dlnorm(xseq, est$mu, est$sigma), ylim = c(0,4e-5), col = "red", type="l", lwd = 2)
      }
      if (input$distr == "expo") {  
        est <- exponential(aps)
        lines(xseq, dexp(xseq, est$lambda), col = "red", type="l",lwd = 2)
      }
      if (input$distr == "gam") {  
        est <- gam(aps, input$method)
        lines(xseq, dgamma(xseq, shape = est$shape, scale=est$scale), col = "red" , type="l", lwd = 2)
      }
      if (input$distr == "p3") {  
        est <- p3(aps, input$method)
        lines(xseq, PearsonDS::dpearsonIII(xseq, scale = est$scale, shape = est$shape, location = est$location), col = "red" , type="l", lwd = 2)
      }
      if (input$distr == "lp3") {  
        est <- lp3(log(aps), input$method)
        if(input$method == "mle") {
          lines(density(rlogpearson(200000,a=est$scale,b=est$shape, c=est$location)), col = "red" , type="l", lwd = 2)
        } else {
          lines(xseq, dlpearsonIII(xseq, skew = est$skew, sdlog = est$sd, meanlog = est$mean), col = "red" , type="l", lwd = 2)
        }
      }
      if (input$distr == "gum") {  
        est <- gumb(aps, input$method)
        lines(xseq, evd::dgumbel(xseq, loc = est$location, scale = est$scale), col = "red" , type="l", lwd = 2)
      }
      if (input$distr == "wei") {  
        est <- weibull(aps, input$method)
        lines(xseq, dweibull3(xseq, shape = est$shape, thres = est$threshold, scale = est$scale), col = "red" , type="l", lwd = 2)
      }
      lines(density(aps), col = "blue", lwd = 2)
    })
    
    output$quantilePlot <- renderPlot({
      return_periods <- (dim(data_sorted)[1]+1)/data_sorted$rank    # Weibull plotting position   
      plot(log(return_periods), data_sorted$aps, xlab = "log return period", ylab = "Discharge (cfs)")
      
      if (input$distr == "norm") {  
        est <- normal(aps, input$method, log = FALSE)
        xt <- qnorm(1-(1/return_periods), est$mu, est$sigma)
        lines(log(return_periods), xt, col=2)
      }
      if (input$distr == "lognorm") {  
        est <- normal(aps, input$method, log = TRUE)
        xt <- qlnorm(1-(1/return_periods), est$mu, est$sigma)
        lines(log(return_periods), xt, col=2)
      }
      if (input$distr == "expo") {  
        est <- exponential(aps)
        xt <- log(return_periods)*(1/est$lambda)
        lines(log(return_periods), xt, col=2)
      }
      if (input$distr == "gam") {  
        est <- gam(aps, input$method)
        z <- qnorm(1-(1/return_periods))
        k <- skew(aps)/6
        kt <- z+(z^2-1)*k + (1/3)*(z^3-6*z)*k^2-(z^2-1)*k^3+z*k^4+(1/3)*k^5
        xt <- est$shape*est$scale + kt*sqrt(est$shape*(est$scale)^2)
        lines(log(return_periods), xt, col=2)
      }
      if (input$distr == "p3") {  
        est <- p3(aps, input$method)
        z <- qnorm(1-(1/return_periods))
        k <- skew(aps)/6
        kt <- z+(z^2-1)*k + (1/3)*(z^3-6*z)*k^2-(z^2-1)*k^3+z*k^4+(1/3)*k^5
        xt <- est$location + est$shape*est$scale + kt*sqrt(est$shape*(est$scale)^2)
        lines(log(return_periods), xt, col=2)
      }
      if (input$distr == "lp3") {  
        est <- p3(log(aps), input$method)
        z <- qnorm(1-(1/return_periods))
        k <- skew(log(aps))/6
        kt <- z+(z^2-1)*k + (1/3)*(z^3-6*z)*k^2-(z^2-1)*k^3+z*k^4+(1/3)*k^5
        xt <- est$location + est$shape*est$scale + kt*sqrt(est$shape*(est$scale)^2)
        lines(log(return_periods), exp(xt), col=2)
      }
      if (input$distr == "gum") {  
        est <- gumb(aps, input$method)
        xt <- est$location - est$scale*log(-log(1-1/return_periods))
        lines(log(return_periods), xt, col=2)
      }
      if (input$distr == "wei") {  
        est <- weibull(aps, input$method)
        xt <- est$threshold + est$scale*(log(return_periods))^(1/est$shape)
        lines(log(return_periods), xt, col=2)
      }
      lines(density(aps), col = "blue", lwd = 2)
    })
    
    # output$info <- renderText({
    #   xy <- function(e) {
    #     if(is.null(e)) return("NULL\n")
    #     paste0("x=", round(e$x, 1), " y=", round(e$y, 1), "\n")
    #   }
    #   paste0(
    #     "hover: ", xy(input$plot_hover)
    #   )
    # })
}

# Run the application 
shinyApp(ui = ui, server = server)

