library(readxl)
library(sqldf)
library(devtools)
library(fitdistrib)
library(EnvStats)
library(nsRFA)
library(MASS)
library(weibullness)
library(FAdist)
library(Lmoments)
library(evd)
library(smwrBase)
library(lmomco)
library(e1071)
# -------------------------------------------------------------------------------------------
# IMPORT & PROCESS DATA
# -------------------------------------------------------------------------------------------
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

# high outlier checking 
lmu <- mean(log(data_processed$aps)); ls <- sd(log(data_processed$aps))
upr <- exp(lmu + 3.049*ls)

# remove outliers
data_cleaned <- data_processed[data_processed$aps < upr,]

# low outlier checking 
lmu <- mean(log(data_cleaned$aps)); ls <- sd(log(data_cleaned$aps))
lwr <- exp(lmu - 3.049*ls)

# remove outliers
data_cleaned <- data_processed[data_processed$aps > lwr,]

# rank observations in decreasing order
data_ranked <- data_cleaned[order(data_cleaned$aps, decreasing = TRUE),]
data_ranked["rank"] <- seq(1,dim(data_cleaned)[1])
data_ranked <- data_ranked[order(data_ranked$Year),]
attach(data_ranked)

# order by rank
data_sorted <- data_ranked[order(data_ranked$rank),]

# for plotting
xseq <- seq(0,100000,10)
# -------------------------------------------------------------------------------------------
# ESTIMATION
# -------------------------------------------------------------------------------------------
# estimation for (log)normal distribution
normal <- function(x, method, log = FALSE) {
  if (log == TRUE) { x <- log(x) }
  if ((method == "mle") || (method == "mme")) {
    list("mu" = mean(x), "sigma" = sqrt(mean((x - mean(x))^2)))
  } else if (method == "pwme") { 
    lmom <- Lmoments(x)
    list("mu" = lmom[1], "sigma" = lmom[2]*sqrt(pi)) }
}
# -------------------------------------------------------------------------------------------
# estimation for  exponential distribution
exponential <- function(x) {
  list("lambda" = 1/mean(x))
}
# -------------------------------------------------------------------------------------------
# estimation for  gamma distribution
gam <- function(x, method) {
  if ((method == "mle") || (method == "mme")) {
    params <-  egamma(x, method = method)$parameter
    # fitdist(x, "gamma", method, start=list(shape=0.5, scale=0.5))
    list("shape" = params[1], "scale" = params[2])
  } else if (method == "pwme") { 
    params <- gamma_Lmom(x)$estimate
    list("shape" = params[1], "scale" = 1/params[2])}
}
# -------------------------------------------------------------------------------------------
# estimation for Pearson 3 distribution
p3 <- function(x, method) {
  if (method == "mle") {
    params <-  pearson_mle(x)$estimate
  } else if (method == "mme") { 
    params <- pearson_mom(x)$estimate
  } else if (method == "pwme") { 
    params <- pearson_Lmom(x)$estimate
  }
    list("location" = params[1], "scale" = params[2], "shape" = params[3])
}
# -------------------------------------------------------------------------------------------
# estimation for Log-Pearson 3 distribution
lp3 <- function(x, method) {
  if (method == "mle") {
    ll <- function(x, par){
    if(par[1]>0 &  par[2]>0 & par[3]> -min(x)) {
      return( -sum(dgamma(x+par[3],shape=par[1],scale=par[2],log=TRUE))) } else {
        return(Inf)
      }
    }
   params <- optim(x = x, c(3,5,5), nll, method="BFGS")$par
   list("location" = params[3], "scale" = params[2], "shape" = params[1])
  } else if (method == "mme") { 
    list("mean" = mean(x), "sd" = sd(x), "skew" = skew(x))
  } else if (method == "pwme") { 
    params <- parpe3(lmoms(x), checklmom=TRUE)$para
    list("mean" = params[1], "sd" = params[2], "skew" = params[3])
  }
}
rlogpearson <- function(n,a,b,c) return( exp(rgamma(n,shape=a,scale=b) - c) )
# -------------------------------------------------------------------------------------------
# estimation for Gumbel distribution
gumb <- function(x, method) {
  params <-  eevd(x, method = method)$parameter
  list("location" = params[1], "scale" = params[2])
}
# -------------------------------------------------------------------------------------------
# estimation for Weibull distribution
weibull <- function(x, method) {
  if (method == "mle") {
    params <-  weibull.mle(x)
    list("shape" = params$shape, "scale" = params$scale, "threshold" = params$threshold)
  } else if (method == "mme"){
    mu1 <- moment(x); mu2 <- moment(x, order = 2); mu3 <- moment(aps, order = 3); mu4 <- moment(aps, order = 2);
    alpha <- (log(2))/(log(mu1-mu2)-log(mu2-mu4))
    mu <- (mu1*mu4 - mu2^2)/(mu1+mu4-2*mu2)
    beta <- (mu1 - mu)/gamma(1+1/alpha)
    # memp <- function(x, order) {mean(x^order)}
    # fitdist(x, "weibull", method="mme", order=c(1,2,3), memp=memp, lower=c(0,0))
    list("shape" = alpha, "scale" = beta, "threshold" = mu)
  } else if (method == "pwme") {
    params <- parwei(lmoms(x),checklmom=TRUE)$para
    list("shape" = params[3], "scale" = params[2], "threshold" = -params[1])
  }
}
