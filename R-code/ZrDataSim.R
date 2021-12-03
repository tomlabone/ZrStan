###############################################################################
#
# Generate bioassay and model data (individual and population) 
# for stable Zirconium and output to CmdStan compatible input files
# using simulated data.
#
# Tom LaBone
# September 2, 2021
#
###############################################################################

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = 1)
library(rethinking)
set.seed(1957)

rm(list=ls(all=TRUE)) 

###############################################################################
#Define biokinetic model
###############################################################################

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

#log of prior covariance matrix -- no correlations
gsd <- 1.5
sigma <- diag(log(gsd)^2,nt)

###############################################################################
#Bioassay data for individuals
###############################################################################

setwd("~/Dropbox/LaBone/GitHub/data")

data.all <- read.csv("ZrSimData-612.csv")
run <- 50
data.in <- data.all[data.all$run==run,]

setwd("~/Dropbox/LaBone/GitHub/subjects")

subjects <- unique(data.in$subject)
N <- length(subjects)

for (i in 1:N) {
  Sub <- subjects[i]
  data.sub <- data.in[data.in$subject == Sub, ]
  M <- data.sub$M
  u <- data.sub$u
  T <- data.sub$T
  dT <- data.sub$dT
  nu <- nrow(data.sub[data.sub$dT != 0,])
  np <- nrow(data.sub[data.sub$dT == 0,])
  n <- np + nu
  
  mu_v <- log(data.sub$vol[i])
  sigma_v <- log(1.1)
  
  #df for robust regression
  df <- 4
  
  #relative uncertainty
  ru <- abs(u/M)
  
  #Output data file for cmdstan
  fname <- gsub(" ","",Sub,fixed = TRUE)
  stan_rdump(
    c(
      "np",
      "nu",
      "n",
      "M",  
      "u",
      "ru",
      "T",    
      "dT",
      "H",
      "nt",
      "nc",
      "theta",
      "sigma",
      "mu_v",
      "sigma_v",
      "df"
    ),
    file = paste(fname,".data.R",sep="")
  )
  
  #Output init files for cmdstan
  V <- 3.4 
  Kt <- theta
  
  stan_rdump(
    c("V","Kt"),
    file = paste(fname,".init.R",sep="")
  )
}

###############################################################################
# Dump simulated population data for cmdstan run
###############################################################################

#read in individual Stan dump files 
dataP <- list()

subFile <- paste(gsub(" ","",subjects[1],fixed = TRUE),".data.R",sep="")
tmp1 <- read_rdump(subFile)
dataM <- tmp1[c("nc","nt","H","df","theta")] 
tmp2 <- tmp1[c("T","dT","M","u","np","nu","n","mu_v","sigma_v")]
dataP[length(dataP)+1] <- list(tmp2)
for (i in 2:length(subjects)) {
  sub <- gsub(" ","",subjects[i],fixed = TRUE)
  subFile <- paste(sub,".data.R",sep="")
  tmp1 <- read_rdump(subFile)
  tmp2 <- tmp1[c("T","dT","M","u","np","nu","n","mu_v","sigma_v")]
  dataP[length(dataP)+1] <- list(tmp2)
}

#These are model parameters -- not person specific
nc <- dataM$nc #the number of compartments in the biokinetic model
nt <- dataM$nt #the number of transfer rate constants in the biokinetic model
H <- dataM$H  #adjacency matrix

#hyper parameters for mean theta of population rate constants
mu_theta <- dataM$theta
gsd <- 1.5
sigma_theta <- diag(log(gsd)^2,nt)

#hyper parameters for tau as log sd
mu_tau <- log(1.55)
sigma_tau <- log(1.1)

#hyper parameter for Omega
Eta <- 10

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Set up data structures for Stan

ntot <- 0
n <- NULL
sub <- NULL
T <- NULL
dT <- NULL
M <- NULL
u <- NULL
mu_v <- NULL
sigma_v <- NULL

# Number of sunjects in population model
N <- 5

for(i in 1:N) {
  n <- c(n,dataP[[i]]$n)               #total number of bioassay per subject
  ntot <- ntot + n[i]                  #total number of bioassay for 16 subjects
  sub <- c(sub,rep(i,dataP[[i]]$n))    #subject array
  T <- c(T,dataP[[i]]$T)               #time bioassay was performed
  dT <- c(dT,dataP[[i]]$dT)            #time interval for urine samples
  M <- c(M,dataP[[i]]$M)               #bioassay measurements
  u <- c(u,dataP[[i]]$u)               #uncertainty in measurements
  mu_v <- c(mu_v,dataP[[i]]$mu_v)     #logmean for blood-volume lognormal prior
  sigma_v <- c(sigma_v,dataP[[i]]$sigma_v) #logsd for blood-volume lognormal prior
}

df <- 4

###############################################################################
# Dump data for cmdstan run
###############################################################################

setwd("~/Dropbox/LaBone/GitHub/cmdstan-home")

stan_rdump(c(
  "ntot",
  "N",
  "n",
  "sub",
  "T",
  "dT",
  "M",
  "u",
  "nt",
  "nc",
  "H",
  "mu_v",
  "sigma_v",
  "mu_theta",
  "sigma_theta",
  "mu_tau",
  "sigma_tau",
  "Eta",
  "df"),
  file = paste("Zr612-",run,".data.R",sep="")
)

#Output init files for cmdstan
V = rep(log(3.4),N)
Kt = matrix(rep(mu_theta,N),nrow=nt,ncol=N)
Omega = diag(nt)
tau = 0.25
theta = mu_theta

stan_rdump(c(
  "V",
  "Kt",
  "Omega",
  "tau",
  "theta"),
  file = paste("Zr612-", run, ".init.R", sep = "")
)

###############################################################################
# Spaghetti plots
###############################################################################

setwd("~/Dropbox/LaBone/GitHub/subjects")

#read in individual Stan dump files 
dataP <- list()

subFile <- paste(gsub(" ","",subjects[1],fixed = TRUE),".data.R",sep="")
tmp1 <- read_rdump(subFile)
dataM <- tmp1[c("nc","nt","H","df","theta")] 
tmp2 <- tmp1[c("T","dT","M","u","np","nu","n","mu_v","sigma_v")]
dataP[length(dataP)+1] <- list(tmp2)
for (i in 2:N) {
  sub <- gsub(" ","",subjects[i],fixed = TRUE)
  subFile <- paste(sub,".data.R",sep="")
  tmp1 <- read_rdump(subFile)
  tmp2 <- tmp1[c("T","dT","M","u","np","nu","n","mu_v","sigma_v")]
  dataP[length(dataP)+1] <- list(tmp2)
}

setwd("~/Dropbox/LaBone/GitHub/plots")

png("PlasmaSpagSim.png",width = 650,height = 480)
plot(
  dataP[[1]]$T[dataP[[1]]$dT==0],
  dataP[[1]]$M[dataP[[1]]$dT==0],
  type="l",
  log="xy",
  ylim=c(0.0001,1),
  xlim=c(0.001,200),
  xlab="Time After Injection (days)",
  ylab="Concentration of Zr in Plasma (fraction of injection/kg plasma)"
)
for(i in 1:N) {
  lines(
    dataP[[i]]$T[dataP[[i]]$dT==0],
    dataP[[i]]$M[dataP[[i]]$dT==0]
  )
  points(
    dataP[[i]]$T[dataP[[i]]$dT==0],
    dataP[[i]]$M[dataP[[i]]$dT==0],
    pch=19,
    cex=0.6
  )
}
dev.off()

png("UrineSpagSim.png",width = 650,height = 480)
plot(
  dataP[[1]]$T[dataP[[1]]$dT>0],
  dataP[[1]]$M[dataP[[1]]$dT>0],
  type="l",
  log="xy",
  ylim=c(0.00001,0.1),
  xlim=c(0.2,200),
  xlab="Time After Injection (days)",
  ylab="Excretion of Zr in Urine (fraction of injection/hour)"
)
for(i in 1:N) {
  lines(
    dataP[[i]]$T[dataP[[i]]$dT>0],
    dataP[[i]]$M[dataP[[i]]$dT>0],
  )
  points(
    dataP[[i]]$T[dataP[[i]]$dT>0],
    dataP[[i]]$M[dataP[[i]]$dT>0],
    pch=19,
    cex=0.5
  )
}
dev.off()






