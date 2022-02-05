setwd("C:/Users/Shiyu Peng/Desktop/QF5203 Project PengShiyu")
library(rugarch)
library(rmgarch)
library(VineCopula)
library(MASS)

# market value of equity
dep = as.matrix(read.csv("r-dep.csv"))[-1,-1]
ins = as.matrix(read.csv("r-ins.csv"))[-1,-1]
bro = as.matrix(read.csv("r-bro.csv"))[-1,-1]
oth = as.matrix(read.csv("r-oth.csv"))[-1,-1]
W = cbind(dep,ins,bro,oth)
W[is.na(W)] <- 0

# log returns of firms
rdep = diff(log(dep))
rins = diff(log(ins))
rbro = diff(log(bro))
roth = diff(log(oth))
r = cbind(rdep,rins,rbro,roth)
r[is.infinite(r)] <- 0
r[is.na(r)] <- 0

# balance-sheet
bal = as.matrix(read.csv(file = "bal.csv"))[,-1]
A = bal[2:53,]
L = bal[55:106,]
OUTS = bal[108:159,]
MV = bal[161:212,]

# market
mar = as.matrix(read.csv(file = "mac.csv"))[-1,2]
mar = diff(log(mar))
mar[is.infinite(mar)] <- 0
mar[is.na(mar)] <- 0

# season number
m = 52

# company number
c = 76

# threshold
cth = -0.1

#individual company srisk
srisk = matrix(0, m,c)
LRMES = matrix(0, m,c)

# trading date in a season
M = length(mar)/m

# day
h = 20

# k 
k = 0.08

# excutive time
time = matrix(0, c)


# srisk-normal
erf <- function(x) {2 * pnorm(x * sqrt(2)) - 1}

for(i in 1:c){

tim = Sys.time()

try(for (j in 3:(m-1)){
  ri = r[1:round(M*j), i]
  rm = mar[1:round(M*j)]
  dat = cbind(ri,rm)
  
  garch = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                     mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                     distribution.model = "norm")
  garchm = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                      distribution.model = "norm")
  
  gar1 =  ugarchfit(spec = garch, ri)
  sigmar = gar1@fit$sigma[round(M*j)]
  gar2 =  ugarchfit(spec = garchm, rm)
  sigmam = gar2@fit$sigma[round(M*j)]
  
  spec = dccspec(uspec = multispec(c(garch, garchm)), VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
                 lag.criterion = c("AIC"), external.regressors = NULL, 
                 robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500), 
                 dccOrder = c(1,1), model = c("DCC"), groups = rep(1, length(uspec@spec)), 
                 distribution = c("mvnorm"), start.pars = list(), fixed.pars = list()) 
  
  dcc = dccfit(spec, data = dat, out.sample = 0, solver = "solnp", solver.control = list(), 
               fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE))
  rho = dcc@mfit$H[1,2,length(dat)/2]

  beta = rho*sigmar/sigmam
  
  LRMES[j,i] = -exp(h*(beta^2*sigmam^2+(1-rho^2)*sigmar^2)/2)*pnorm((beta*log(1+cth)-h*beta^2*sigmam^2)/sqrt(h)/beta/sigmam)/(0.5+0.5*erf(beta*log(1+cth)/sqrt(2)/sqrt(h)/beta/sigmam))+1
  MVOE = as.numeric(OUTS[j,i])*as.numeric(MV[j,i])/10^6
  LVG = (as.numeric(L[j,i])+MVOE)/MVOE
  srisk[j,i] = MVOE*(k*LVG+(1-k)*LRMES[j,i]-1)
  print(srisk[j,i])
  },silent = TRUE)
time[i] = Sys.time() - tim
print(time[i])
}


write.csv(LRMES, file = "lrmes-nor.csv")
write.csv(time, file = "srisk-nor-time.csv")
write.csv(srisk, file = "srisk-nor.csv")
