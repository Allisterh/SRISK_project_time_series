setwd("C:/Users/Shiyu Peng/Desktop/QF5203 Project PengShiyu")


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

# k 
k = 0.08


## overall
srisk = as.matrix(read.csv("srisk-cop.csv"))[-1,-1]
srisk = as.matrix(read.csv("srisk-nor.csv"))[-1,-1]
srisk = as.matrix(read.csv("srisk-cop(normal).csv"))[-1,-1]
srisk = as.matrix(read.csv("srisk-gar.csv"))[-1,-1]
srisk[is.na(srisk)]<-0

srisk1 = matrix(0, m)
srisk2 = matrix(0, m)
srisk3 = matrix(0, m)
srisk4 = matrix(0, m)
srisk_all = matrix(0, m)

for(k in 1:(m-1)){
  srisk1[k] = sum(max(srisk[k,1:ncol(dep)],0))
  srisk2[k] = sum(max(srisk[k,(ncol(dep)+1):(ncol(dep)+ncol(ins))],0))
  srisk3[k] = sum(max(srisk[k,(ncol(dep)+ncol(ins)+1):(ncol(dep)+ncol(ins)+ncol(bro))],0))
  srisk4[k] = sum(max(srisk[k,(ncol(dep)+ncol(ins)+ncol(bro)+1):(ncol(dep)+ncol(ins)+ncol(bro)+ncol(oth))],0))
}

for(k in 1:(m-1)){
  srisk_all[k] <- srisk1[k]+srisk2[k]+srisk3[k]+srisk4[k]
}
srisk_all<-srisk_all[1:(length(srisk_all)-2)]

x = unique(2000+(1:50)*0.25)
plot(x,srisk_all,type = "l",main = "srisk (Dynamic Rotated Gumbel Copula)",xlab = "year")
plot(x,srisk_all,type = "l",main = "srisk (Static Normal)",xlab = "year")
plot(x,srisk_all,type = "l",main = "srisk (Dynamic Normal Copula)",xlab = "year")

##ranking
srisk_rank = matrix(0, 1,c)
for(k in 1:c){
  srisk_rank[k] = sum(max(srisk[1:(m-1),k],0))
}
write.csv(srisk_rank,file = "srisk-rank-cop.csv")
write.csv(srisk_rank,file = "srisk-rank-nor.csv")
write.csv(srisk_rank,file = "srisk-rank-cop(normal).csv")
rank_cop=as.matrix(read.csv("srisk-rank-cop.csv"))
rank_nor=as.matrix(read.csv("srisk-rank-nor.csv"))
rank_cop_normal=as.matrix(read.csv("srisk-rank-cop(normal).csv"))
order(rank_cop,decreasing=T)
order(rank_nor,decreasing=T)
order(rank_cop_normal,decreasing=T)

##timely
library(lmtest)
sriskcop = matrix(0, m)
srisknor = matrix(0, m)
sriskcopnor = matrix(0, m)
sriskgar = matrix(0, m)
sriskcop = diff(log(srisk_all))
srisknor = diff(log(srisk_all))
sriskcopnor = diff(log(srisk_all))
sriskgar = diff(log(srisk_all))

sriskcop[is.infinite(sriskcop)]<-0
srisknor[is.infinite(srisknor)]<-0
sriskcopnor[is.infinite(sriskcopnor)]<-0
sriskcop[is.na(sriskcop)]<-0
srisknor[is.na(srisknor)]<-0
sriskcopnor[is.na(sriskcopnor)]<-0

grangertest(sriskcop ~ sriskcopnor, order = 1)

##plot
sriskcop = matrix(0, m-2)
srisknor = matrix(0, m-2)
sriskcopnor = matrix(0, m-2)
sriskgar = matrix(0, m-2)
sriskcop = srisk_all
srisknor = srisk_all
sriskcopnor = srisk_all
sriskgar = srisk_all

x = unique(2000+(1:50)*0.25)
plot(x,sriskcop,type = "l",xlab = "year",col="green",ylab = "Srisk",ylim=c(0,310000))
par(new=TRUE)
plot(x,srisknor,axes = FALSE,type = "l",col="red",xlab = "",ylab = "",ylim=c(0,310000))
par(new=TRUE)
plot(x,sriskcopnor,type = "l",col="blue",xlab = "",ylab = "",ylim=c(0,310000))
par(new=TRUE)
plot(x,sriskgar,type = "l",col="blue",xlab = "",ylab = "",ylim=c(0,310000))
par(new=TRUE)