setwd("~/Desktop/Maths Project/code")
# setwd("~/Year 4/Project/code")

library(readxl)
library(devtools)
install_github("NVE/fitdistrib")
# install_github("NVE/FlomKart")
# library("FlomKart")
library(fitdistrib)
library(sqldf)
library(EnvStats)
library(MASS)
library(weibullness)
library(FAdist)
library(fitdistrplus)
library(Lmoments)
library(evd)

# read in data
# data <- read_excel('data/Pearl River - USGS02486000.xlsx', sheet = 1)
data <- read_excel('data/USGS02486100.xlsx', sheet = 1)
attach(data)

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
upr <- mu - 3.049*s

# remove outlier (more than 3 sds above mean)
data_processed <- data_processed[aps < upr,]
attach(data_processed)
plot(Year, aps, type="l")
hist(aps, breaks=15, main = "Distribution of annual peak discharge", xlab="Annual peak streamflow (cf/s)")

# rank observations in decreasing order
data_ranked <- data_processed[order(aps, decreasing = TRUE),]
data_ranked["rank"] <- seq(1,dim(data_processed)[1])
data_ranked <- data_ranked[order(Year),]
attach(data_ranked)

xseq <- seq(0,100000,10)

# -------------------------------------------------------------------------------------------
# estimation for (log)normal distribution
normal <- function(x, method, log = FALSE) {
  if (log == TRUE) { x <- log(x) }
  if ((method == "mle") || (method == "mme")) {
    list("mu" = mean(x), "sigma" = sqrt(mean((x - mean(x))^2)))
  } else if (method == "pwme") { 
    lmom <- Lmoments(x,rmax=2)
    list("mu" = lmom[1], "sigma" = lmom[2]*sqrt(pi)) }
}

# Normal - MLE
norm_mle <- normal(aps, "mle", log = FALSE)
plot(xseq, dnorm(xseq, norm_mle$mu, norm_mle$sigma),ylim = c(0,4e-5), main="Normal (MLE)", col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Normal - MOM
norm_mom <- normal(aps, "mme", log = FALSE)
plot(xseq, dnorm(xseq, norm_mom$mu, norm_mom$sigma), ylim = c(0,4e-5), main="Normal (MOM)", col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Normal - PWM
norm_pwm <- normal(aps, "pwme", log = FALSE)
plot(xseq, dnorm(xseq, norm_pwm$mu, norm_pwm$sigma), ylim = c(0,4e-5), main="Normal (PWM)", col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Log-normal - MLE
lognorm_mle <- normal(aps, "mle", log = TRUE)
plot(xseq[xseq>0], dlnorm(xseq[xseq>0], lognorm_mle$mu, lognorm_mle$sigma), main="Log-normal (MLE)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Log-normal - MOM
lognorm_mom <- normal(aps, "mme", log = TRUE)
plot(xseq[xseq>0], dlnorm(xseq[xseq>0], lognorm_mom$mu, lognorm_mom$sigma), main="Log-normal (MOM)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Log-normal - PWM
lognorm_pwm <- normal(aps, "pwme", log = TRUE)
plot(xseq[xseq>0], dlnorm(xseq[xseq>0], lognorm_pwm$mu, lognorm_pwm$sigma), main="Log-normal (PWM)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# -------------------------------------------------------------------------------------------
# estimation for  exponential distribution
exponential <- function(x, method) {
    list("lambda" = 1/mean(x))
}
exp_mle <- exponential(aps, "mle")
exp_mom <- exponential(aps, "mme")
plot(xseq, dexp(xseq, exp_mle$lambda), main="Exponential (MLE and MOM)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# -------------------------------------------------------------------------------------------
# estimation for  gamma distribution
gam <- function(x, method) {
  if ((method == "mle") || (method == "mme")) {
    params <-  egamma(x, method = method)$parameter
    list("shape" = params[1], "scale" = params[2])
  } else if (method == "pwme") { 
    params <- gamma_Lmom(x)$estimate
    list("shape" = params[1], "scale" = 1/params[2])}
}
# Gamma - MLE
gamma_mle <- gamma(aps, "mle")
plot(xseq, dgamma(xseq, shape = gamma_mle$shape, scale=gamma_mle$scale), main="Gamma (MLE)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Gamma - MOM
gamma_mom <- gamma(aps, "mme")
plot(xseq, dgamma(xseq, shape = gamma_mom$shape, scale=gamma_mom$scale), main="Gamma (MOM)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# -------------------------------------------------------------------------------------------
# estimation for gumbel distribution
gumbel <- function(x, method) {
  params <-  eevd(x, method = method)$parameter
  list("location" = params[1], "scale" = params[2])
}

# Gumbel - MLE
gumbel_mle <- gumbel(aps, "mle")
plot(xseq, dgumbel(xseq, loc = as.numeric(gumbel_mle$location), scale=as.numeric(gumbel_mle$scale)), main="Gumbel (MLE)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Gumbel - MOM
gumbel_mom <- gumbel(aps, "mme")
plot(xseq, dgumbel(xseq, location = gumbel_mom$location, scale=gumbel_mom$scale), main="Gumbel (MOM)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# Gumbel - PWM
gumbel_pwm <- gumbel(aps, "pwme")
plot(xseq, dgumbel(xseq, location = gumbel_pwm$location, scale=gumbel_pwm$scale), main="Gumbel (PWM)",  ylim = c(0,4e-5), col = 1, , type="l", xlab = "", ylab = "density")
lines(density(aps), col = 2)

# -------------------------------------------------------------------------------------------
# estimation for weibull distribution NOT VERIFIED
weibull <- function(x, method) {
  if (method == "mle") {
    params <-  weibull.mle(x)
    list("shape" = params$shape, "scale" = params$scale, "threshold" = params$threshold)
  } else if (method == "mme"){
      print("not available")
  } else if (method == "pwme") {
      print("not available") 
    }
}
weibull_mle <- weibull(aps, "mle")

# Weibull - MLE
plot(xseq, dweibull3(xseq,shape = weibull_mle$shape, scale = weibull_mle$scale, thres = weibull_mle$threshold), ylim = c(0, 4e-5), col = 2, type='l',main = "Weibull (MLE)", xlab = "", ylab = "density")
lines(density(aps))
#________________________________________________________________
# pearson III MLE
nll <- function(x, p) {
  alpha <- p[1]
  beta <- p[2]
  epsilon <- p[3]
  n <- length(x)
  if(alpha>0 & epsilon <= min(x)) {
  # n*log(alpha) + n*log(base::gamma(beta)) - (beta-1)* sum(log(x-epsilon)) + n*(beta-1)*log(alpha)+sum((x-epsilon)/alpha)
    -sum(dgamma(x+epsilon, scale = alpha, shape = beta, log=TRUE))
  } else { 
    return(Inf)
  }
}
params_p3 <- optim(x = aps, p = c(0.5,0.5,0.5), fn = nll, method = "BFGS")

plot(xseq, dpearsonIII(xseq, scale = params_p3[1], shape = params_p3[2], location = params_p3[3]), col = 2, type='l', main="Pearson 3 (MLE)", ylab = "density", xlab="")
lines(density(aps))

# log-pearson III MLE
rlogpearson <- function(n,a,b,c) return( exp(rgamma(n,shape=a,scale=b) - c) )
ll <- function(x, par){
 if(par[1]>0 &  par[2]>0 & par[3]> -min(x)) {
    return( -sum(dgamma(x+par[3],shape=par[1],scale=par[2],log=TRUE))) } else {
      return(Inf)
    }
}
params_lp3 <- optim(x = log(aps), c(3,5,5),nll, method="BFGS")$par

plot(density(rlogpearson(100000,b=params_lp3[1],a=params_lp3[2], c=params_lp3[3])), col = 3, ylim = c(0, 4e-5), xlim = c(-1e4, 1e05), main="Log-Pearson 3 (MLE)", ylab = "density", xlab="")
lines(density(aps))

plot(density(rlogpearsonIII(100000, a = est$shape, b = est$scale, c = est$location)), col = "red" , type="l", lwd = 2))

params_lp3 <- pearson_Lmom(log(aps))$estimate
plot(xseq[xseq>min(log(aps))], dlgamma3(xseq[xseq>min(log(aps))], scale=params_lp3[1],shape=params_lp3[2], thres=params_lp3[3]), col = 3, ylim = c(0, 4e-5), xlim = c(-1e4, 1e05), main="Log-Pearson 3 (MLE)", ylab = "density", xlab="")
lines(density(aps))

#________________________________________________________
# plots
data_sorted <- data_ranked[order(data_ranked$rank),]
return_periods <- (dim(data_sorted)[1]+1)/data_sorted$rank 
probs <- 1/return_periods
plot(log(return_periods), data_sorted$aps)

xt_norm <- function(return_period, params) { params[1] + qnorm(1-(1/return_period), 0, 1)*params[2] }
xt_lnorm <- function(return_period, params) { exp(params[1] + qnorm(1-(1/return_period), 0, 1)*params[2]) }
lines(log(return_periods),xt_norm(return_periods, params = c(mean(aps), sd(aps))), col=2)
lines(log(return_periods),xt_lnorm(return_periods, ), col=3)

model <- lm(aps[1:100] ~ qnorm(1-(1/return_period[1:100])))
model$coefficients <- c(mean(aps), sd(aps))
predictions <- predict(model, data.frame(qnorm(1-(1/seq(1,100,1)))), interval = "confidence")
lines(log(seq(5,101,0.25)), x_t(seq(5,101,0.25))+qt(0.975,384)*sqrt((1+0.5*qnorm(1-(1/seq(5,101,0.25))))*var(aps)/382))


# gof
ad.test(aps, null = "pgamma", shape=3.50, scale = 9083)
-2*sum(dnorm(aps, mean=mean(aps), sd=sd(aps), log=TRUE))+2*2
-2*sum(dnorm(aps, mean=mean(aps), sd=sd(aps), log=TRUE))+2*log(n)
