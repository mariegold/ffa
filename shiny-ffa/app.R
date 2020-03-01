library(shiny)
library(readxl)

# read in data
data <- read_excel('data/Pearl River - USGS02486000.xlsx', sheet = 1)
attach(data)

# sort by date and convert to dataframe
data_processed <- data.frame(data[order(DecYear),])
colnames(data_processed)[2] <- "aps"

# remove period prior to 1900 (not enough data)
data_processed <- data_processed[data_processed$Year > 1899,]

# drop minimum value in years with multiple values
data_processed <- sqldf("select Year, max(aps) as aps from data_processed group by Year")

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Flood Frequency Analysis"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("bins",
                     "Number of bins:",
                     min = 1,
                     max = 50,
                     value = 30)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
      x    <- faithful[, 2] 
      bins <- seq(min(x), max(x), length.out = input$bins + 1)
      
      # draw the histogram with the specified number of bins
      plot(data_processed$Year, data_processed$aps, type="l")
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

