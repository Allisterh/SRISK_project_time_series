setwd("C:/Users/Shiyu Peng/Desktop/QF5203 Project PengShiyu")
library(rugarch)
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

# monte-c
mc = 500

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

# back test
btest = matrix(0, length(mar),c)
btestm = matrix(0, length(mar),c)

# srisk-Copula
for (i in 1:c){
  tim = Sys.time()
try(for (j in 3:(m-1)){
    ri = r[1:round(M*j), i]
    rm = mar[1:round(M*j)]
    dat = cbind(ri,rm)
    garch = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                       mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                       distribution.model = "std")
    garchm = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
                        mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                        distribution.model = "std")
    gar1 =  ugarchfit(spec = garch, ri)
    sigmar = gar1@fit$sigma
    gar2 =  ugarchfit(spec = garchm, rm)
    sigmam = gar2@fit$sigma
    # dynamic copula
    stdri = ri/sigmar
    stdrm = rm/sigmam
    FitT = fitdistr(stdri, 't')
    dfri = FitT$estimate[3]
    FitT = fitdistr(stdrm, 't')   # marginal dist
    dfrm = FitT$estimate[3]
    delta = matrix(0, round(M*j)+h)
    abs = matrix(0, round(M*j)+h)
    aveabs = matrix(0, round(M*j)+h)
    for(t in (round(M*j)-50):round(M*j)){
      dcri = r[1:t,i]
      dcrm = mar[1:t]
      n = length(dcri)
      uri = rank(dcri)/(n+1)
      urm = rank(dcrm)/(n+1)
      try<-try(BiCopEst(uri,urm,family=14),silent=TRUE)
      if('try-error' %in% class(try)){delta[t] = 0 }
      else{Bi = BiCopEst(uri,urm,family=14)
      delta[t] = Bi[2]}
      abs[t] = abs(pt(ri[t]/sigmar[t],df = dfri) - pt(rm[t]/sigmam[t],df = dfrm))
      if(t>(round(M*j)-41-1)){aveabs[t] = mean(abs[(t-9):t])}
    }
    delta = as.numeric(delta)
    deltasub = delta[(round(M*j)-40):round(M*j)] - 1
    deltalag = delta[(round(M*j)-41):(round(M*j)-1)]
    ave = aveabs[(round(M*j)-40):round(M*j)]
    reg = lm(deltasub ~ (1+ave+deltalag)^2)
    for(t in 1:h){
      abs[round(M*j)+t] = abs(pt(rt(1,df = dfri),df = dfri) - pt(rt(1,df = dfrm),df = dfrm)) # uniform
      aveabs[round(M*j)+t] = mean(abs[(round(M*j)+t-9):(round(M*j)+t)])
    }
    for(t in (round(M*j)+1):(round(M*j)+h))
    {delta[t] = 1+(reg$coefficients[1]+reg$coefficients[2]*aveabs[t]+reg$coefficients[3]*delta[t-1])^2
    if(delta[t]>100){delta[t]=100}}
    simr = matrix(0, h, mc)
    simm = matrix(0, h, mc)
    fore = ugarchforecast(gar1, n.ahead = h)
    fsigmar = fore@forecast$sigmaFor
    fore = ugarchforecast(gar2, n.ahead = h)
    fsigmam = fore@forecast$sigmaFor
    for(x in 1:h){
      sim = BiCopSim(mc, family = 14, par = delta[round(M*j)+x])
      for(t in 1:mc){
        simr[x,t] = fsigmar[x]*qt(sim[t,1], df = dfri)
        simm[x,t] = fsigmam[x]*qt(sim[t,2], df = dfrm)
      }
      btest[round(M*j)+x,i] = quantile(simr[x,], 0.05)
      btestm[round(M*j)+x,i] = quantile(simm[x,], 0.05)
    }
    ind = matrix(0, mc,1)
    R1h = matrix(0, mc,1)
    for (t in 1:mc) {ind[t] = ( (exp(sum(simm[,t])) - 1) <= cth); R1h[t] = exp(sum(simr[,t])) - 1}
    LRMES[j,i] = -sum(R1h*ind)/sum(ind)
    MVOE = as.numeric(OUTS[j,i])*as.numeric(MV[j,i])/10^6
    LVG = (as.numeric(L[j,i])+MVOE)/MVOE
    srisk[j,i] = MVOE*(k*LVG+(1-k)*LRMES[j,i]-1)
    print(srisk[j,i])
  },silent = TRUE)
  time[i] = Sys.time() - tim
  print(time[i])
}

write.csv(LRMES, file = "lrmes-cop.csv")
write.csv(btest, file = "srisk-cop-btest.csv")
write.csv(btestm, file = "srisk-cop-btestm.csv")
write.csv(time, file = "srisk-cop-time.csv")
write.csv(srisk, file = "srisk-cop.csv")
