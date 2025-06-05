####REAL DATA ANALYSIS####
##Data Extraction
source("MPGIG_INGARCH_FN (git).R")
datt = read.csv("Cannabis_example.csv")
datt = as.matrix(datt[,-1])

ts.plot(datt,col=1:2)
acf(datt)
pacf(datt)

####Implementation of EM algorithm
#old_INGARCH_EM(dat = data,init=initial value if you have, model=list(past_mean=,past_obs=),trace=T,method="BFGS")
em.1.12.old = old_INGARCH_EM(datt,init=NULL,model=list(past_mean=NULL,past_obs=c(1,12)),trace=T,method="BFGS")
em.1.12.poi = new_INGARCH_Poi(datt,init=NULL,model=list(past_mean=NULL,past_obs=c(1,12)),method="BFGS")

#Diagnostic
acf(resid(em.1.12.old))
acf(fitted(em.1.12.old))

##PIT histograms
J=10
tmp.old = pit_mSichel_old(datt,em.1.12.old,1) # For the first series
tmp.old2 = pit_mSichel_old(datt,em.1.12.old,2) # For the second series
pit_1.old = pit_mar(c(1:J)/J,tmp.old) ;pit_1.old = c(pit_1.old[1],diff(pit_1.old))*10
pit_2.old = pit_mar(c(1:J)/J,tmp.old2);pit_2.old = c(pit_2.old[1],diff(pit_2.old))*10
tmp.poi = Poi_u(datt,em.1.12.poi,1)
tmp.poi2 = Poi_u(datt,em.1.12.poi,2)
pit_poi_1 = pit_mar(c(1:J)/J,tmp.poi) ;pit_poi_1 = c(pit_poi_1[1],diff(pit_poi_1))*10
pit_poi_2 = pit_mar(c(1:J)/J,tmp.poi2);pit_poi_2 = c(pit_poi_2[1],diff(pit_poi_2))*10

par(mfrow=c(2,2))
hist(rep(seq(0.05,0.95,by=0.1),round(pit_1.old*1000)),breaks = 10,
     main="MPGIG (MNC)",freq=F,ylab="",xlab="",ylim=c(0,2.5),cex.main=2)
abline(h=1,lwd=2,lty=2,col="blue")
hist(rep(seq(0.05,0.95,by=0.1),round(pit_2.old*1000)),breaks = 10,
     main="MPGIG (GNC)",freq=F,ylab="",xlab="",ylim=c(0,2.5),cex.main=2)
abline(h=1,lwd=2,lty=2,col="blue")

hist(rep(seq(0.05,0.95,by=0.1),round(pit_poi_1*1000)),breaks = 10,
     main="Poisson (MNC)",freq=F,ylab="",xlab="",ylim=c(0,2.5),cex.main=2)
abline(h=1,lwd=2,lty=2,col="blue")
hist(rep(seq(0.05,0.95,by=0.1),round(pit_poi_2*1000)),breaks = 10,
     main="Poisson (GNC)",freq=F,ylab="",xlab="",ylim=c(0,2.5),cex.main=2)
abline(h=1,lwd=2,lty=2,col="blue")
par(mfrow=c(1,1))





##Parametric Boostrap
set.seed(77)
em.1.12.old.pb = old_ingarch_pb(em.1.12.old,pb_n=500)
em.1.12.old.pb.names= c(expression(phi),expression(alpha),expression(d[1]),expression(d[2]),
                        expression(b[11(1)]),expression(b[21(1)]),expression(b[12(1)]),expression(b[22(1)]),
                        expression(b[11(12)]),expression(b[21(12)]),expression(b[12(12)]),expression(b[22(12)]))
x11();par(mfrow=c(3,4))
for (i in 1:4) hist(em.1.12.old.pb$pb.est[,i],xlab="",ylab="",freq=F,breaks=20,main=em.1.12.old.pb.names[i],cex.main=2.2)
for (i in c(5,7,9,11,6,8,10,12)) hist(em.1.12.old.pb$pb.est[,i],xlab="",ylab="",freq=F,breaks=20,main=em.1.12.old.pb.names[i],cex.main=2.2)

