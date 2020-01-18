setwd("~/Desktop/Maths Project/Jackson, MS")
library(readxl)
library(devtools)
install_github("NVE/fitdistrib")
install_github("NVE/FlomKart")
library("FlomKart")

# read in data
data <- read_excel('Pearl River - USGS02486000.xlsx', sheet = 1)
attach(data)

# sort by date and convert to dataframe
data <- data.frame(data[order(DecYear),])

plot(DecYear, `Annual Peak Streamflow (cfs)`)

# MLE estimation
library(FAmle)
mle(x=data[, 2], dist='weibull',start=c(.1,.1))
