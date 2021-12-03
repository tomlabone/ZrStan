###############################################################################
#
# Posterior Predictive Distribution for Zr individual model using
# canned priors the simulated population model.
#
# Tom LaBone
# September 28, 2021
#
###############################################################################

library(scales)     #transparent points in base plots
library(magicaxis)
library(expm)
library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
library(posterior)
set.seed(1957)

rm(list=ls(all=TRUE))

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
# Functions
###############################################################################

#Calculate total removal rate constants in matrix K
trrc <- function(k, lambda) {
  K <- matrix(0,nrow(k), nrow(k))
  K = k
  for (i in 1:nrow(K)) {
    Ksum = 0
    for (j in 1:nrow(K)) {
      Ksum = Ksum + k[i, j]
    }
    K[i, i] = -Ksum - lambda
    
  }
  return(K)
}

#calculate compartment content using the matrix exponential method
q_me <- function(comp, x, q0, k) {
  a <- matrix(0,nrow(k), nrow(k))
  q <- numeric(length(x))
  for (i in 1:length(x)) {
    a <- expm(k * x[i])
    for (j in 1:nrow(k)) {
      q[i] = q[i] + a[j, comp] * q0[j]
    }
  }
  return(q)
}

q_es <- function(comp, x, q0, k) {
  es <- eigen(t(k))
  K <- es$values
  V <- es$vectors
  M <- solve(V,q0)
  q <- numeric(length(x))
  for(i in 1:length(x)) {
    for (j in 1:dim(k)[1]) {  
      q[i] <- q[i] + M[j]*V[comp,j]*exp(x[i]*K[j])
    }
  }
  return(q)
} 

#calculate content at times x and rate constants kt
q_c <- function(comp, x, kt, H, nc, nt) {
  q0 <- numeric(nc)
  k <- matrix(0,nc, nc)
  content <- numeric(length(x))
  q0[1] <- 1.0
  for (i in 1:nt) {
    k[H[i, 1], H[i, 2]] <- kt[i]
  }
  k <- trrc(k, 0)
  content <- q_me(comp, x, q0, k)
  return(content)
}


pltCI <- function(x,y,lcl,ucl,lcol="red",type="y") {
  if(type=="y") {
    for(i in 1: length(x)) {
      lines(c(x[i],x[i]),c(lcl[i],ucl[i]),col=lcol)
    }
  } else {
    for(i in 1: length(x)) {
      lines(c(lcl[i],ucl[i]),c(y[i],y[i]),col=lcol)
    }
    
  }
}

###############################################################################
#Read in bioassay data and biokinetic model for subject
###############################################################################

setwd("~/Dropbox/LaBone/GitHub/subjects")
sub <- "Zr_10"
run <- 50

source(paste(sub,".data.R",sep=""))

setwd("~/Dropbox/LaBone/GitHub/data")
#load the canned priors
#theta_sigma
#tau_p
#theta_mu
#v_p
source("Zrcan-50.R")

###############################################################################
# Run Stan model
###############################################################################

mu_theta <- theta_mu
sigma_theta <- theta_sigma
mu_tau <- tau_p[1]
sigma_tau <- tau_p[2]
mu_beta <- 0
sigma_beta <- log(1.2)
Eta <- 10
df <- df
beta <- 0
ru <- abs(u/M)

#prior for blood plasma volume is based on weight/sex algorithm and
#not the canned priors
mu_v <- log((14.5*82 + 2035) / 1000)
sigma_v <- log(1.2)

setwd("~/Dropbox/LaBone/GitHub/cmdstan-home/Zr")

stan_rdump(
  c(
    "np",
    "nu",
    "n",
    "T",
    "dT",
    "M",
    "u",
    "ru",
    "nt",
    "nc",
    "H",
    "mu_v",
    "sigma_v",
    "mu_theta",
    "sigma_theta",
    "mu_tau",
    "sigma_tau",
    "mu_beta",
    "sigma_beta",
    "Eta",
    "df"
  ),
  file = paste(sub,".data.R",sep="")
)

v <- mu_v
Kt <- mu_theta
Omega = diag(nt)
tau = tau_p[1]
theta = mu_theta

stan_rdump(
  c(
    "Kt",
    "Omega",
    "tau",
    "theta",
    "beta",
    "v"
  ),
  file = paste(sub,".init.R",sep="")
)

################################################################################
################################################################################
# Regression
###############################################################################
################################################################################

setwd("~/Dropbox/LaBone/GitHub/data")

damod <- read.csv("Zr-ICRP.csv")

#the number of compartments in the model
nc <- max(c(damod$ToNum,damod$FromNum))

#the number of transfer rate constants in the model
nt <- length(damod[,1])

#adjacency matrix
H <- as.matrix(damod[,c(3,5)])

#log of priors for rate constants
theta <- log(as.numeric(damod$K))

#ICRP 89 plasma volume
V <- 3.0

kt <- exp(theta)

m_p <- ( q_c(1,T[1:np],kt,H,nc,nt) + q_c(2,T[1:np],kt,H,nc,nt) ) / V
m_u <- q_c(15,T[(np+1):n],kt,H,nc,nt) - q_c(15,T[(np+1):n]-dT[(np+1):n],kt,H,nc,nt);
m <- c(m_p,m_u)

fit.lm <- lm(M ~ m - 1, weights=(u/M)^2)

summary(fit.lm)
beta <- signif(summary(fit.lm)$coefficients[1],4)
stderr <- signif(summary(fit.lm)$coefficients[2],4)

setwd("~/Dropbox/LaBone/GitHub/plots")

png(paste(sub,"-RegIntake.png",sep=""),width=650,height=480)
plot(
  m,M,
  xlab="Reference Bioassay Function",
  ylab="Zr in Plasma (conc) and Urine (content)",
  pch=16,
  cex=0.7,
  main=sub
)
abline(a=0,b=1,lty=6,col="pink")
lines(m,fitted(fit.lm))
lines(m,predict(fit.lm,interval="confidence")[,1],col="azure4")
lines(m,predict(fit.lm,interval="confidence")[,2],col="azure4")
lines(m,predict(fit.lm,interval="confidence")[,3],col="azure4")
pltCI(m,M,M-1.96*u,M+1.96*u) 
text(grconvertX(0.01,"npc"),grconvertY(0.8,"npc"),
     bquote(beta==.(beta)),pos=4)
text(grconvertX(0.01,"npc"),grconvertY(0.7,"npc"),
     bquote(sigma==.(stderr)),pos=4)
dev.off()

################################################################################
################################################################################
# Run zircon-ind.stan and analyze output
###############################################################################
################################################################################

# Get the run times from the script output file

setwd("~/Dropbox/LaBone/GitHub/cmdstan-home")

ZrScript <- readLines("Zr-ind-10.out")
run.time <- NULL
for (i in 1:length(ZrScript)) {
  if (substr(ZrScript[i],1,2)=="Zr") {
    run.time <- rbind(run.time,substr(ZrScript[i],1,9))
  }
  if (grepl("Total",ZrScript[i],fixed=TRUE)) {
    run.time <- rbind(run.time,substr(ZrScript[i],1,38))
  }
}
run.time <- gsub(" ","",run.time[,1])
run.time <- round(as.numeric(gsub("([0-9]+).*$", "\\1",run.time))/3600/24,3)
run.time <- paste(run.time,"days")
run.time

###############################################################################
# get input data and output data for each chain
###############################################################################

setwd("~/Dropbox/LaBone/GitHub/cmdstan-home/Zr")

data.out <- read.csv(
  "Zr_10.csv",
  as.is = TRUE,
  header = TRUE,
  comment.char = "#"
)

fit <- rstan::read_stan_csv("Zr_10.csv")
data.stat <- summarize_draws(fit)

###############################################################################
# Plots

setwd("~/Dropbox/LaBone/GitHub/plots")

png(paste(sub,"-RhatESS.png",sep=""),width = 650,height = 480)
plot(
  data.stat$rhat,data.stat$ess_bulk,
  ylab="Bulk ESS",
  xlab="R-hat",
  log="xy",
  col = alpha("black",1/5)
)
dev.off()

Beta <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "Be"])
Beta.mu <- signif(data.stat$mean[substr(data.stat$variable, 1, 2) == "Be"],4)
dens <- density(Beta)
db <- approxfun(dens$x, dens$y)
x <- rnorm(100000,beta,stderr)
denx <- density(x)
dx <- approxfun(denx$x, denx$y)

png(paste(sub,"-Intake.png",sep=""),width=650,height=480)
plot(
  dens,
  xlab="Beta",
  main=sub,
  ylim=c(0,max(max(dens$y),max(denx$y)))
)
lines(
  c(quantile(Beta,probs=0.025),quantile(Beta,probs=0.025)),
  c(0,db(quantile(Beta,probs=0.025))),
  col="black",
  lty=6
)
lines(
  c(quantile(Beta,probs=0.975),quantile(Beta,probs=0.975)),
  c(0,db(quantile(Beta,probs=0.975))),
  col="black",
  lty=6
)
lines(denx,col="red")
lines(
  c(quantile(x,probs=0.975),quantile(x,probs=0.975)),
  c(0,dx(quantile(x,probs=0.975))),
  col="red",
  lty=6
)
lines(
  c(quantile(x,probs=0.025),quantile(x,probs=0.025)),
  c(0,dx(quantile(x,probs=0.025))),
  col="red",
  lty=6
)
abline(v=1,lty=6,col="azure4")
dev.off()

#get the actual blood plasma volume
setwd("~/Dropbox/LaBone/GitHub/data")
tmp <- read.csv("ZrSimData-612.csv")
tmp <- tmp[tmp$run==run & tmp$subject==sub,]
V.true <- tmp$vol[1]
setwd("~/Dropbox/LaBone/GitHub/plots")

V <- as.matrix(data.out[, substr(names(data.out), 1, 1) == "V"])
dens <- density(V)
dv <- approxfun(dens$x, dens$y)
png(paste(sub,"-Vol.png",sep=""),width=650,height=480)
plot(
  density(V),
  xlab="Plasma Volume (L)",
  main=sub
)
abline(v=V.true,col="red",lty=6)
lines(
  c(quantile(V,probs=0.025),quantile(V,probs=0.025)),
  c(0,dv(quantile(V,probs=0.025))),
  col="black",
  lty=6
)
lines(
  c(quantile(V,probs=0.975),quantile(V,probs=0.975)),
  c(0,dv(quantile(V,probs=0.975))),
  col="black",
  lty=6
)
dev.off()

m <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "m."])
M.hat <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "M_"])
M.hat.mu <- apply(M.hat,2,mean)
uclM <- apply(M.hat,2,quantile,probs=0.975)
lclM <- apply(M.hat,2,quantile,probs=0.025)
uclm <- apply(m,2,quantile,probs=0.975)
lclm <- apply(m,2,quantile,probs=0.025)
m.mu <- apply(m,2,mean)

png(paste(sub,"-0.png",sep=""),width=650,height=480)
plot(m.mu,M,
     main=sub,
     xlab="Intake Retention Fraction",
     ylab="Zr in Plasma (conc) and Urine (content)",
     log="",
     cex=0.7,pch=16
)
lines(m.mu,M.hat.mu)
lines(m.mu,m.mu*quantile(Beta,probs=0.975),col="azure4")
lines(m.mu,m.mu*quantile(Beta,probs=0.025),col="azure4")
text(grconvertX(0.01,"npc"),grconvertY(0.8,"npc"),
     bquote(beta==.(Beta.mu)),pos=4)
abline(a=0,b=1,col="pink")
pltCI(m.mu,M,M-1.96*u,M+1.96*u,"black") 
dev.off()

png(paste(sub,"-1.png",sep=""),width=650,height=480)
plot(m.mu,M,
     pch = 19,
     cex = 0.5,
     col = "red",
     xlab="Reference Bioassay Function",
     ylab="Zr in Plasma (conc) and Urine (content)",
     main=paste(sub," Posterior Predictive")
)
lines(m.mu, M.hat.mu)
pltCI(m.mu, M.hat.mu, lclM, uclM)
pltCI(m.mu, M.hat.mu, lclm, uclm,type="x")
points(m.mu,M,pch = 19,cex = 0.5,col = "black")
dev.off()

png(paste(sub,"-2.png",sep=""),width=650,height=480)
plot(m.mu,M,
     pch = 19,
     cex = 0.7,
     col = "red",
     xlim=c(0,0.4),
     ylim=c(0,0.4),
     xlab="Reference Bioassay Function",
     ylab="Zr in Plasma (conc) and Urine (content)",
     main=paste(sub," Posterior Predictive")
)
for (j in 1:length(m[, 1])) {
  points(m[j, ],M.hat[j, ],pch = 19,cex = 0.5,col = alpha("azure4",1/20))
}
lines(m.mu, M.hat.mu)
#pltCI(m.mu,M,M-1.96*u,M+1.96*u,"black") 
pltCI(m.mu, M.hat.mu, lclM, uclM,lcol="red")
pltCI(m.mu, M.hat.mu, lclm, uclm,type="x",lcol="red")
points(m.mu,M,pch = 19,cex = 0.7,col = "black")
dev.off()

################################################################################
# K and Tau Parameter plots
################################################################################

tmp <- data.out
Kt <- (tmp[,substr(colnames(tmp), 1, 2) == "Kt"])
Kt_mu <- apply(Kt,2,mean)

Kt_mode <- numeric(nt)
for(i in 1:nt) {
  dens <- density(Kt[,i],n=2^10)
  index <- which(dens$y==max(dens$y))
  Kt_mode[i] <- dens$x[index]
}
Kt_sd <- apply(Kt,2,sd)

y.lim <- c(0,3)
y.mode <- 0.6
y <- 2.5
png(paste(sub,"-K.png",sep=""),width=650,height=480)
plot(
  density(Kt[,1],n=2^10),
  xlim=c(-5,3),
  ylim=y.lim,
  xlab="Log Rate Constant",
  col="azure4",
  main=paste(sub,"Rate Matrix K Marginal Posteriors")
)
for(i in 2:length(Kt)) {
  lines(density(Kt[,i],n=2^10),col="azure4")
}
doi <- NULL
for(i in 1:27) {
  dens <- density(Kt[,i],n=2^10)
  if(max(dens$y)>y.mode) {
    doi <- c(doi,i)
    lines(density(Kt[,i],n=2^10),col="red")
  }
}
abline(v=Kt_mode[doi],lty=6,col="pink")
x <- Kt_mode[doi]
text(x,y,
     K.names[doi],
     cex=0.7,srt=90
)
dev.off()

K.cor <- cor(Kt)
colnames(K.cor) <- K.names2
row.names(K.cor) <- K.names2
png(paste("CorrPlot-",sub,".png"),width = 650,height = 480)
corrplot::corrplot(K.cor,tl.cex=0.7)
dev.off()

K.cov <- cov(Kt)
colnames(K.cov) <- K.names2
row.names(K.cov) <- K.names2
png(paste("CovPlot-",sub,".png"),width = 650,height = 480)
corrplot::corrplot(K.cov,is.corr=FALSE,tl.cex=0.7)
dev.off()

png(paste(sub,"-Kgsd.png",sep=""),width=650,height=480)
plot(
  exp(Kt_sd),
  xlab="Parameter Index",
  ylab="Geometric Standard Deviation",
  ylim=c(1,3),
  type="n",
)
points(1:27,exp(Kt_sd),pch=16)
text(1:27,exp(Kt_sd)+0.2,label=K.names,cex=0.7,srt=90)
dev.off()

tau.pos <- as.matrix(data.out[, substr(names(data.out), 1, 2) == "ta"])
png(paste(sub,"-Tau.png",sep=""),width=650,height=480)
plot(
  density(exp(tau.pos)),
  xlab="Tau",
  main=paste(sub,"Tau")
)
text(grconvertX(0.5,"npc"),grconvertY(0.9,"npc"),
     paste("gm =",signif(exp(mean(tau.pos)),4)),pos=4)
text(grconvertX(0.5,"npc"),grconvertY(0.85,"npc"),
     paste("gsd =",signif(exp(sd(tau.pos)),4)),pos=4)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Individual plots of K

pdf(paste("Kdensity-",sub,".pdf",sep=""))

for(j in 1:length(K.names)) {
  plot(
    density(Kt[,j]),
    xlim=c(-11,3),
    ylim=c(0,3.5),
    xlab="Log Rate Constant",
    main="Rate Matrix K Marginal Posteriors"
  )
  for(i in 1:length(K.names)) {
    lines(density(Kt[,i]),col="azure4")
  }
  lines(density(Kt[,j]),col="red",lwd=2)
  abline(v=Kt_mu[j],lty=6,col="pink")
  text(Kt_mu[j],2.5,
       K.names[j],
       cex=1
  )
  text(
    grconvertX(0.01,"npc"),grconvertY(0.9,"npc"),
    paste("Rhat =",signif(data.stat$rhat[j],7)),pos=4,col="red"
  )
  text(
    grconvertX(0.01,"npc"),grconvertY(0.95,"npc"),
    paste("n_eff =",round(data.stat$ess_bulk[j],0)),pos=4,col="red"
  )
  text(
    grconvertX(0.7,"npc"),grconvertY(0.95,"npc"),
    paste("mean =",round(data.stat$mean[j],7)),pos=4,col="red"
  )
  text(
    grconvertX(0.7,"npc"),grconvertY(0.9,"npc"),
    paste("sd =",round(data.stat$sd[j],7)),pos=4,col="red"
  )
  text(
    grconvertX(0.4,"npc"),grconvertY(0.95,"npc"),
    paste("gm =",round(exp(data.stat$mean[j]),7)),pos=4,col="red"
  )
  text(
    grconvertX(0.4,"npc"),grconvertY(0.9,"npc"),
    paste("gsd =",round(exp(data.stat$sd[j]),7)),pos=4,col="red"
  )
}
dev.off()

################################################################################
#MAP estimates
################################################################################

setwd("~/Dropbox/LaBone/GitHub/cmdstan-home/Zr")

data.opt <- read.csv(
  paste("Zr_10opt.csv",sep=""),
  as.is = TRUE,
  header = TRUE,
  comment.char = "#"
)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("~/Dropbox/LaBone/GitHub/plots")

Beta <- signif(as.numeric(data.opt[, substr(names(data.opt), 1, 2) == "Be"]),4)
V <- signif(as.numeric(data.opt[, substr(names(data.opt), 1, 1) == "V"]),4)
m <- as.numeric(data.opt[, substr(names(data.opt), 1, 2) == "m."])
M.hat <- as.numeric(data.opt[, substr(names(data.opt), 1, 1) == "M"])

png(paste(sub,"-MAPIntake0.png",sep=""),width=650,height=480)
plot(
  m,M,
  xlab="Reference Bioassay Function",
  ylab="Zr in Plasma (conc) and Urine (content)",
  pch=16,
  cex=0.7,
  main=sub
)
abline(a=0,b=1,lty=6,col="pink")
lines(m,m*Beta)
pltCI(m,M,M-1.96*u,M+1.96*u,"black") 
text(grconvertX(0.01,"npc"),grconvertY(0.8,"npc"),
     bquote(beta==.(Beta)),pos=4)
text(grconvertX(0.01,"npc"),grconvertY(0.75,"npc"),
     paste("v =",V,"L"),pos=4)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use MAP estimates of parameters to calculate reference bioassay functions

V <- data.opt[, substr(names(data.opt), 1, 1) == "V"]
kt <- exp(as.numeric(data.opt[, substr(names(data.opt), 1, 2) == "Kt"]))

m_p <- ( q_c(1,T[1:np],kt,H,nc,nt) + q_c(2,T[1:np],kt,H,nc,nt) ) / V
m_u <- q_c(15,T[(np+1):n],kt,H,nc,nt) - q_c(15,T[(np+1):n]-dT[(np+1):n],kt,H,nc,nt);
m <- c(m_p,m_u)

fit.map <- lm(M ~ m - 1,weights=(u/M)^2)
summary(fit.map)

beta <- signif(summary(fit.map)$coefficients[1],4)
stderr <- signif(summary(fit.map)$coefficients[2],4)

png(paste(sub,"-MAPIntake1.png",sep=""),width=650,height=480)
plot(
  m,M,
  xlab="Reference Bioassay Function",
  ylab="Zr in Plasma (conc) and Urine (content)",
  pch=16,
  cex=0.7,
  main=sub
)
abline(a=0,b=1,lty=6,col="pink")
lines(m,fitted(fit.map))
lines(m,predict(fit.map,interval="confidence")[,1],col="azure4")
lines(m,predict(fit.map,interval="confidence")[,2],col="azure4")
lines(m,predict(fit.map,interval="confidence")[,3],col="azure4")
pltCI(m,M,M-1.96*u,M+1.96*u,"black") 
text(grconvertX(0.01,"npc"),grconvertY(0.8,"npc"),
     bquote(beta==.(beta)),pos=4)
text(grconvertX(0.01,"npc"),grconvertY(0.7,"npc"),
     bquote(sigma==.(stderr)),pos=4)
dev.off()

