setwd("~/Desktop/Maths Project/code")
# setwd("~/Year 4/Project/code")

library(readxl)
library(devtools)
install_github("NVE/fitdistrib")
install_github("NVE/FlomKart")
library("FlomKart")
library(sqldf)
library(fitdistrib)

# read in data
data <- read_excel('data/Pearl River - USGS02486000.xlsx', sheet = 1)

# sort by date and convert to dataframe
data_processed <- data.frame(data[order(DecYear),])
colnames(data_processed)[2] <- "aps"

# remove period prior to 1900 (not enough data)
data_processed <- data_processed[data_processed$Year > 1899,]

# drop minimum value in years with multiple values
data_processed <- sqldf("select Year, max(aps) as aps from data_processed group by Year")

attach(data_processed)

# visualisation
plot(Year, aps, type="l")

# outlier checking
mu <- mean(aps); s <- sd(aps)
upr <- mu + 3*s

# remove outlier (more than 3 sds above mean)
data_processed <- data_processed[aps < upr,]
attach(data_processed)
plot(Year, aps, type="l")

# rank observations in decreasing order
data_ranked <- data_processed[order(aps, decreasing = TRUE),]
data_ranked["rank"] <- seq(1,dim(data_processed)[1])
data_ranked <- data_ranked[order(Year),]
attach(data_ranked)

# MLE estimation ??
library(FAmle)
plot(mle(x=aps, dist='weibull',start=c(.1,.1)))

library(MASS)
# fitdistr(aps, "gamma") DOESN'T WORK
