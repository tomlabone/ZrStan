###############################################################################
#
# Analysis of Zirconium Bioassay Data with Stan - posterior predictive plots
# Uses simulated data. 
#
# Tom LaBone
# August 15, 2021
#
###############################################################################

rm(list=ls(all=TRUE)) 

library(MASS)
library(rstan)
library(scales) #alpha
library(mclust)
library(chron)
library(posterior) 

pltCI <- function(x,y,lcl,ucl) {
  for(i in 1: length(x)) {
    lines(c(x[i],x[i]),c(lcl[i],ucl[i]),col="red")
  }
}
 
K.names <- c(
  "K[bld1,bld2]",
  "K[bld1,liv0]",
  "K[bld1,kid]",
  "K[bld1,STO]",
  "K[bld1,ST1]",
  "K[bld1,ubc]",
  "K[bld1,SI]",
  "K[bld1,ts]",
  "K[bld1,cs]",
  "K[bld2,bld1]",
  "K[liv0,SI]",
  "K[liv0,bld1]",
  "K[liv0,liv1]",
  "K[liv1,bld1]",
  "K[kid,bld1]",
  "K[STO,bld1]",
  "K[ST1,bld1]",
  "K[ts,bld1]",
  "K[ts,tv]",
  "K[tv,bld1]",
  "K[cs,bld1]",
  "K[cs,cv]",
  "K[cv,urn]",
  "K[ubc,urn]",
  "K[SI,col]",
  "K[col,fec]",
  "K[SI,bld1]"
)

K.names2 <- c(
  "bld1,bld2",
  "bld1,liv0",
  "bld1,kid",
  "bld1,ST0",
  "bld1,ST1",
  "bld1,ubc",
  "bld1,SI",
  "bld1,TS",
  "bld1,CS",
  "bld2,bld1",
  "liv0,SI",
  "liv0,bld1",
  "liv0,liv1",
  "liv1,bld1",
  "kid,bld1",
  "STO,bld1",
  "ST1,bld1",
  "TS,bld1",
  "TS,TV",
  "TV,bld1",
  "CS,bld1",
  "CS,CV",
  "CV,bld1",
  "ubc,urn",
  "SI,col",
  "col,fec",
  "SI,bld1"
)

###############################################################################
# Get the run times from the script output file

setwd("~/Dropbox/LaBone/Github/cmdstan-home")

run <- 50
ZrScript <- readLines("Zr-20210612-50.out")
run.time <- NULL
for (i in 1:length(ZrScript)) {
  if (substr(ZrScript[i],1,2)=="Zr") {
    run.time <- rbind(run.time,substr(ZrScript[i],1,9))
  }
  if (grepl("Total",ZrScript[i],fixed=TRUE)) {
    run.time <- rbind(run.time,substr(ZrScript[i],1,38))
  }
}
run.time <- strsplit(run.time[,1]," ")
run.time <- run.time[[1]][16]
run.time <- round(as.numeric(run.time)/3600/24,3)
run.time <- paste(run.time,"days")
run.time

###############################################################################
# get input data and output data for each chain
###############################################################################

source("Zr20210612-50.data.R")
M.all <- M
u.all <- u 

setwd("~/Dropbox/LaBone/Github/cmdstan-home/Zr")

data1 <- read.csv(
  "Zr-20210612-50.csv",
  header=TRUE,
  comment.char = "#"
)
fit <- rstan::read_stan_csv("Zr-20210612-50.csv")

stats <- summarize_draws(fit)

data.out <- data1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#multiple chains are combined and then analyzed like this

#fit <- sflist2stanfit(list(fit1,fit2,fit3,fit4,fit5,fit6))
#stats <- summarize_draws(fit)
#data.out <- rbind(data1,data2,data3,data4,data5,data6)

###############################################################################
# get subject names

setwd("~/Dropbox/LaBone/Github/data")

data.in <- read.csv("ZrSimData-20210612.csv")
subjects <- sort(unique(data.in$subject))
subjects

#We  only used 5 subjects for this example
subjects <- subjects[1:5]

###############################################################################
# set up datasets for analysis

V.all <- as.matrix(data.out[substr(names(data.out), 1, 2) == "V."])
v.all <- exp(V.all)
V.mu <- apply(V.all,2,mean)
V.sd <- apply(V.all,2,sd)
v.mu <- exp(V.mu)
v.sd <- exp(V.sd)

m.all <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "m."])
M.hat.all <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "M_"])

theta <- data.frame(as.matrix(data.out[, substr(names(data.out), 1, 2) == "th"]))
tau <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "ta"])
Tau <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "Ta"])

###############################################################################
# posterior predictive plots

setwd("~/Dropbox/LaBone/Github/plots")

for (i in 1:length(subjects)) {
  Sub <- gsub(" ","",subjects[i],fixed=TRUE)
  m <- m.all[,which(sub==i)]
  M.hat <- M.hat.all[,which(sub==i)]
  v <- v.all[,1]
  M <- M.all[which(sub==i)] 
  u <- u.all[which(sub==i)]
  
  v.mu <- mean(v)
  M.hat.mu <- apply(M.hat,2,mean)
  ucl <- apply(M.hat,2,quantile,probs=0.975)
  lcl <- apply(M.hat,2,quantile,probs=0.025)
  
  ucl.gum <- M + 2*u
  lcl.gum <- M - 2*u
  
  m.mu <- apply(m,2,mean)
  
  png(paste(Sub, "-Pop.png", sep = ""),width = 650,height = 480)
  plot(m.mu,M,
       pch = 19,
       cex = 0.5,
       col = "red",
       ylim = c(-0.1, 0.5),
       xlim = c(0, 0.40),
       xlab = "Reference Bioassay Function",
       ylab = "Observed and Predicted Bioassay",
       main=Sub
  )
  for (j in 1:length(m[, 1])) {
    points(m[j, ],M.hat[j, ],pch = 19,cex = 0.5,col = alpha("azure4",1/7))
  }
  lines(m.mu, M.hat.mu)
  pltCI(m.mu, M.hat.mu, lcl, ucl)
  
  points(m.mu,M,pch = 19,cex = 0.7,col = "black")
  for(i in 1: length(m.mu)) {
    lines(c(m.mu[i],m.mu[i]),c(M[i]-2*u[i],M[i] + 2*u[i]))
  }
  dev.off()
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#empirical blood plasma

png("Plasma-Pop.png",width = 650,height = 480)
plot(
  density(V.all[,1]),
  xlim=c(0.75,2.5),
  ylim=c(0,20),
  main="",
  xlab="Log of Blood Plasma Volume"
)
for(i in 1:dim(V.all)[2]) {
  lines(density(V.all[,i]))
}
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#empirical correlation and covariance of theta

theta.cor <- cor(theta)
colnames(theta.cor) <- K.names2
row.names(theta.cor) <- K.names2
png("CorrPlot-Pop.png",width = 650,height = 480)
corrplot::corrplot(theta.cor,tl.cex=0.7)
dev.off()

theta.cov <- cov(theta)
colnames(theta.cov) <- K.names2
row.names(theta.cov) <- K.names2
png("CovPlot-Pop.png",width = 650,height = 480)
corrplot::corrplot(theta.cov,is.corr=FALSE,tl.cex=0.7)
dev.off()

###############################################################################
# MCMC Daignostics
###############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Trace plots for parameters tracked by Stan

setwd("~/Dropbox/LaBone/Github/diag")

### Theta
tmp <- data.out[,substr(names(data.out),1,2) == "th"]
stat <- stats[substr(stats$variable,1,2) == "th", ]
for (i in 1:nt) {
  cname <- names(tmp)[i]
  png(paste(cname, "-Pop.png", sep = ""),width = 650,height = 480)
  plot(
    tmp[, i],
    type = "l",
    col = alpha("black", 1 / 2),
    ylab = cname,
    main = K.names[i],
    xlab = "Iteration"
  )
  se <- stat$sd[i] / sqrt(stat$ess_bulk[i])
  rhat <- stat$rhat[i]
  neff <- stat$ess_bulk[i]
  text(
    grconvertX(0.01, "npc"),
    grconvertY(0.85, "npc"),
    paste("SE =", signif(se, 7)),
    pos = 4,
    col = "red"
  )
  text(
    grconvertX(0.01, "npc"),
    grconvertY(0.9, "npc"),
    paste("Rhat =", signif(rhat, 7)),
    pos = 4,
    col = "red"
  )
  text(
    grconvertX(0.01, "npc"),
    grconvertY(0.95, "npc"),
    paste("n_eff =", round(neff, 0)),
    pos = 4,
    col = "red"
  )
  abline(h = stat$mean[i], col = "red")
  print(i)
  dev.off()
}

### Blood plasma volume
tmp <- data.out[,substr(names(data.out),1,1) == "V"]
stat <- stats[substr(stats$variable,1,1) == "V", ]
for (i in 1:N) {
  cname <- names(tmp)[i]
  png(paste(cname, "-Pop.png", sep = ""),width = 650,height = 480)
  plot(
    tmp[, i],
    type = "l",
    col = alpha("black", 1 / 2),
    ylab = cname,
    main = cname[i],
    xlab = "Iteration"
  )
  se <- stat$sd[i] / sqrt(stat$ess_bulk[i])
  rhat <- stat$rhat[i]
  neff <- stat$ess_bulk[i]
  text(
    grconvertX(0.01, "npc"),
    grconvertY(0.85, "npc"),
    paste("SE =", signif(se, 7)),
    pos = 4,
    col = "red"
  )
  text(
    grconvertX(0.01, "npc"),
    grconvertY(0.9, "npc"),
    paste("Rhat =", signif(rhat, 7)),
    pos = 4,
    col = "red"
  )
  text(
    grconvertX(0.01, "npc"),
    grconvertY(0.95, "npc"),
    paste("n_eff =", round(neff, 0)),
    pos = 4,
    col = "red"
  )
  abline(h = stat$mean[i], col = "red")
  print(i)
  dev.off()
}

### tau

stats[substr(stats$variable,1,2) == "Ta", ]
stats[substr(stats$variable,1,2) == "ta", ]

tmp <- data.out[,substr(names(data.out),1,2) == "ta"]
stat <- stats[substr(stats$variable,1,2) == "ta", ]

png("tau-Pop.png",width = 650,height = 480)
plot(
  tmp,
  type="l",
  col = alpha("black",1/2),
  ylab="tau",
  main="tau",
  xlab="Iteration"
)
se <- stat$sd / sqrt(stat$ess_bulk)
rhat <- stat$rhat
neff <- stat$ess_bulk
text(
  grconvertX(0.01, "npc"),
  grconvertY(0.85, "npc"),
  paste("SE =", signif(se, 7)),
  pos = 4,
  col = "red"
)
text(
  grconvertX(0.01, "npc"),
  grconvertY(0.9, "npc"),
  paste("Rhat =", signif(rhat, 7)),
  pos = 4,
  col = "red"
)
text(
  grconvertX(0.01, "npc"),
  grconvertY(0.95, "npc"),
  paste("n_eff =", round(neff, 0)),
  pos = 4,
  col = "red"
)
abline(h = stat$mean, col = "red")
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# R-hat and N.eff

png("ESS-Pop.png",width = 650,height = 480)
plot(
  stats$ess_bulk,stats$ess_tail,
  log="xy",
  col = alpha("black",1/2),
  xlab="Bulk Effective Sample Size",
  ylab="Tail Effective Sample Size"
)
dev.off()

png("RhatESS-Pop.png",width = 650,height = 480)
plot(
  stats$rhat,stats$ess_bulk,
  ylab="Bulk ESS",
  xlab="R-hat",
  log="xy"
)
dev.off()

###############################################################################
# look at the parameter estimates
###############################################################################

setwd("~/Dropbox/LaBone/Github/plots")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#tau

stat <- (stats[substr(stats$variable, 1, 2) == "ta",])
stderr <- stat[4] / sqrt(stat[9])

tau <- exp(as.matrix(data.out[, substr(names(data.out), 1, 2) == "ta"]))
dens <- density(tau)
dt <- approxfun(dens$x, dens$y)
png("Tau-Pop.png",width=650,height=480)
plot(
  density(tau),
  xlab="Tau",
  main=""
)
lines(
  c(quantile(tau,probs=0.025),quantile(tau,probs=0.025)),
  c(0,dt(quantile(tau,probs=0.025))),
  col="black",
  lty=6
)
lines(
  c(quantile(tau,probs=0.975),quantile(tau,probs=0.975)),
  c(0,dt(quantile(tau,probs=0.975))),
  col="black",
  lty=6
)
lines(
  c(quantile(tau,probs=0.5),quantile(tau,probs=0.5)),
  c(0,dt(quantile(tau,probs=0.5))),
  col="black",
  lty=6
)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Theta marginal posteriors

stat <- stats[substr(stats$variable, 1, 2) == "th",]

png("ThetaMarginal-Pop0.png",width = 650,height = 480)
plot(
  density(theta[,1]),
  xlim=c(-11,3),
  ylim=c(0,3),
  xlab="Log Theta",
  main=""
)
for(i in 2:nrow(stat)) {
  lines(density(theta[,i]),col="azure4")
}
for(i in c(1,4,6,8,9,10,16)) {
  lines(density(theta[,i]),col="red")
}
abline(v=stat$mean[c(1,4,6,8,9,10,16)],lty=6,col="pink")
x <- c(stat$mean[c(1,4,6,8,9,10,16)])

y <- c(2.3,2.6)
text(x,y,
     K.names[c(1,4,6,8,9,10,16)],
     cex=0.7,srt=90
)
dev.off()

png("ThetaMarginal-Pop1.png",width = 650,height = 480)
plot(
  density(theta[,1]),
  xlim=c(-4,3),
  ylim=c(0,3),
  xlab="Log Theta",
  main=""
)
for(i in 2:nrow(stat)) {
  lines(density(theta[,i]),col="azure4")
}
for(i in c(1,4,6,8,9,10,16)) {
  lines(density(theta[,i]),col="red")
}
abline(v=stat$mean[c(1,4,6,8,9,10,16)],lty=6,col="pink")
x <- c(stat$mean[c(1,4,6,8,9,10,16)])

y <- c(2.3,2.6)
text(x,y,
     K.names[c(1,4,6,8,9,10,16)],
     cex=0.7,srt=90
)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Individual plots of theta

pdf(paste("Theta-Pop.pdf",sep=""))
for(j in 1:length(K.names)) {
  plot(
    density(theta[,j]),
    xlim=c(-11,3),
    ylim=c(0,3.5),
    xlab="Log Rate Constant",
    main="Theta Marginal Posteriors"
  )
  for(i in 1:length(K.names)) {
    lines(density(theta[,i]),col="azure4")
  }
  lines(density(theta[,j]),col="red",lwd=2)
  abline(v=stat$mean[j],lty=6,col="pink")
  text(stat$mean[j],2.5,
       K.names[j],
       cex=1
  )
  text(
    grconvertX(0.01,"npc"),grconvertY(0.9,"npc"),
    paste("Rhat =",signif(stat$rhat[j],7)),pos=4,col="red"
  )
  text(
    grconvertX(0.01,"npc"),grconvertY(0.95,"npc"),
    paste("n_eff =",round(stat$ess_bulk[j],0)),pos=4,col="red"
  )
  text(
    grconvertX(0.7,"npc"),grconvertY(0.95,"npc"),
    paste("mean =",round(stat$mean[j],7)),pos=4,col="red"
  )
  text(
    grconvertX(0.7,"npc"),grconvertY(0.9,"npc"),
    paste("sd =",round(stat$sd[j],7)),pos=4,col="red"
  )
  text(
    grconvertX(0.4,"npc"),grconvertY(0.95,"npc"),
    paste("gm =",round(exp(stat$mean[j]),7)),pos=4,col="red"
  )
  text(
    grconvertX(0.4,"npc"),grconvertY(0.9,"npc"),
    paste("gsd =",round(exp(stats$sd[j]),7)),pos=4,col="red"
  )
}
dev.off()

png("GSD-theta.png",width=650,height=480)
plot(
  exp(stat$sd),
  xlab="Parameter Index",
  ylab="Geometric Standard Deviation",
  ylim=c(1,2),
  type="n",
)
points(1:27,exp(stat$sd),pch=16)
text(1:27,exp(stat$sd)+0.1,label=K.names,cex=0.7,srt=90)
dev.off()

###############################################################################
# Generate canned priors
###############################################################################

setwd("~/Dropbox/LaBone/Github/data")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Blood plasma volume

V_p <- NULL
for(i in 1:N) {
  fit1 <- fitdistr(V.all[,i],"normal")
  V_p <- rbind(V_p,as.numeric(fit1$estimate))
}
v_p <- c(mean(V_p[,1]),max(V_p[,2]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# theta

fit1 <- mvn("Ellipsoidal", data = theta)
fit1$loglik

#parameter estimates
theta_mu <- as.numeric(fit1$parameters$mean)
theta_mu

theta_sigma <- matrix(
  fit1$parameters$variance$sigma,
  nrow = length(theta_mu),
  ncol = length(theta_mu)
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# tau

tmp <- data.out[,substr(names(data.out),1,2) == "ta"]
stat <- stats[substr(stats$variable,1,2) == "ta", ]
fit1 <- fitdistr(tmp,"normal")
tau_p <-as.numeric(fit1$estimate)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# dump canned priors

stan_rdump(
  c(
    "v_p",
    "tau_p",
    "theta_mu",
    "theta_sigma"
  ),
  file = paste("Zrcan-",run,".R",sep="")
)

