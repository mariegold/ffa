source('fitting.R')

Jackson = "USGS02486000"
Edinburg = "USGS02482000"
Carthage = "USGS02482550"
Lena = "USGS02483500"
Rockport = "USGS02488000"
Monticello = "USGS02488500"
Columbia = "USGS02489000"
Bogalusa = "USGS02489500"

data <- read_excel(paste0("data/", Jackson, ".xls"), sheet = 1) 
out <- preprocessing(data)
data_processed <- out$processed
data_sorted <- out$sorted
aps <- out$aps
return_periods <- (dim(data_sorted)[1]+1)/data_sorted$rank

fit <- function(distr, method) {
  if (distr == "norm") {
    estim <- normal(aps, method, log = FALSE, rp = return_periods, plotx = xseq)
  }
  if (distr == "lognorm") {
    estim <- normal(log(aps), method, log = TRUE, rp = return_periods, plotx = xseq)
  }
  if (distr == "expo") {
    estim <- exponential(aps, rp = return_periods, plotx = xseq)
  }
  if (distr == "gam") {
    estim <- gam(aps, method, rp = return_periods, plotx = xseq)
  }
  if (distr == "p3") {
    estim <- p3(aps, method, rp = return_periods, plotx = xseq)
  }
  if (distr == "lp3") {
    estim <- p3(log(aps), method, rp = return_periods, log=TRUE)
    estim$xt <- exp(estim$xt)
    if (method == "mle") {
      estim$likelihood <- lp3(log(aps), method)$likelihood
    }
  }
  if (distr == "gum") {
    estim <- gumb(aps, method, rp = return_periods, plotx = xseq)
  }
  if (distr == "wei") {
    estim <- weibull(aps, method, rp = return_periods, plotx = xseq)
  }
  estim
}

fit("lnorm", "pwme")$par


# fitdistr(aps, "gamma", start=list(scale = 10, shape=1.5))$sd


plot(log(return_periods), data_sorted$aps, xlab = "Log return period", ylab = "Discharge (cfs)", 
     main = "Estimated peak discharge vs log return period")
estim <- fit("norm", "pwme")
lines(log(return_periods), estim$xt, col=2, lwd=2)
estim <- fit("lognorm", "pwme")
lines(log(return_periods), estim$xt, col=3, lwd=2)
estim <- fit("expo", "mle")
lines(log(return_periods), estim$xt, col="orange", lwd=2)
estim <- fit("gam", "mle")
lines(log(return_periods), estim$xt, col=5, lwd=2)
estim <- fit("p3", "mle")
lines(log(return_periods), estim$xt, col=6, lwd=2)
estim <- fit("lp3", "mle")
lines(log(return_periods), estim$xt, col=4, lwd=2)
estim <- fit("gum", "mme")
lines(log(return_periods), estim$xt, col=8, lwd=2)
estim <- fit("wei", "pwme")
lines(log(return_periods), estim$xt, col=9, lwd=2)
legend("bottomright", col = c(2,3,"orange",5,6,4,8,9), lty = rep(1,8), c("Normal","Log-normal", "Exponential", "Gamma", "Pearson III", "Log-Pearson III", "Gumbel", "Weibull"))





