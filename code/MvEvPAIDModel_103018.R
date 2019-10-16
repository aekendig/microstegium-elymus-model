# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)

# set working directory
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/REU/Data")

# import data (all data combined up to post harvest)
dat=read.csv("MvEv_GerminationInfection_PostHarvest.csv",header=T)
head(dat)

# subset data
unique(dat$Date)
unique(dat$Treatment)
# dates selected based on stabilization diagrams from REU analysis
mvDat=subset(dat,Date%in%c(20180721,20180722,20180723,20180724)&Treatment=="None"&SpPresent=="Mv")
evDat=subset(dat,Date%in%c(20180721,20180722,20180723,20180724)&Treatment=="None"&SpPresent=="Ev")

# germination rates
mvDat$germ.mv=with(mvDat,GermMv/50)
evDat$germ.ev=with(evDat,GermEv/50)
(germ.mv=mean(mvDat$germ.mv))
(germ.ev=mean(evDat$germ.ev))

# seed production rates
(lambda.mv=2500/10*12.5/20) # Wilson et al. 2015
lambda.ev=435 # Stevens 1957

# initial numbers
N0.mv=1
N0.ev=1

# intraspecific competition
alpha.mm=seq(0,2,length=10)
alpha.ee=seq(0,2,length=10)

# interspecific competition
alpha.me=seq(0,2,length=10)
alpha.em=seq(0,2,length=10)

# years
simtime=100

# Ev by itself
Ni.ev=rep(NA,simtime)
Ni.ev[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev[t]=Ni.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[8]*germ.ev*Ni.ev[t-1])
}

# plot
df.ev=data.frame(N=Ni.ev,time=1:simtime)
ggplot(df.ev,aes(x=time,y=N))+geom_point()

# Mv invades Ev
N.mv=rep(NA,simtime)
N.ev=rep(NA,simtime)

N.mv[1]=N0.mv
N.ev[1]=Ni.ev[simtime]

# introduce Mv population
for(t in 2:simtime){
	N.ev[t]=N.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[8]*germ.ev*N.ev[t-1]+alpha.em[10]*germ.mv*N.mv[t-1])
	N.mv[t]=N.mv[t-1]*germ.mv*lambda.mv/(1+alpha.mm[2]*germ.mv*N.mv[t-1]+alpha.me[3]*germ.ev*N.ev[t-1])
}

# plot
df=data.frame(N=c(N.ev,N.mv),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime))
ggplot(df,aes(x=time,y=N))+geom_point(aes(colour=species))

# infection reduces seed production of Mv
lambdaI.mv=lambda.mv*0.26 # Flory et al. 2011

# Infected Mv population
NI.mv=rep(NA,simtime)
NI.mv[1]=N0.mv

# Healthy Ev population
NH.ev=rep(NA,simtime)
NH.ev[1]=Ni.ev[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NH.ev[t]=NH.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[8]*germ.ev*NH.ev[t-1]+alpha.em[10]*germ.mv*NI.mv[t-1])
	NI.mv[t]=NI.mv[t-1]*germ.mv*lambdaI.mv/(1+alpha.mm[2]*germ.mv*NI.mv[t-1]+alpha.me[3]*germ.ev*NH.ev[t-1])
}

# plot
dfI=data.frame(N=c(NH.ev,NI.mv),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),status="Infected")
df$status="Healthy"
df2=merge(df,dfI,all=T)
ggplot(df2,aes(x=time,y=N))+geom_point(aes(colour=species))+facet_wrap(~status)

## Add-ons
# declining lambda for infected status (accumulates disease and fecundity declines)
# parameter combination of alphas that leads to native persistence vs. decline
# seed reduction value that leads to native persistence vs. decline

# declining lambda
spreadTime=5
arriveTime=10
lambdadI.mv=c(rep(lambda.mv,arriveTime),seq(lambda.mv,lambdaI.mv,length=spreadTime),rep(lambdaI.mv,simtime-(spreadTime+arriveTime)))

# Infected Mv population
NdI.mv=rep(NA,simtime)
NdI.mv[1]=N0.mv

# Healthy Ev population
NdH.ev=rep(NA,simtime)
NdH.ev[1]=Ni.ev[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev[t]=NdH.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[5]*germ.ev*NdH.ev[t-1]+alpha.em[10]*germ.mv*NdI.mv[t-1])
	NdI.mv[t]=NdI.mv[t-1]*germ.mv*lambdadI.mv[t]/(1+alpha.mm[5]*germ.mv*NdI.mv[t-1]+alpha.me[6]*germ.ev*NdH.ev[t-1])
	NdH.ev[t]=ifelse(NdH.ev[t]<1,0,NdH.ev[t])
	NdI.mv[t]=ifelse(NdI.mv[t]<1,0,NdI.mv[t])
}

# plot
dfdI=data.frame(N=c(NdH.ev,NdI.mv),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),status="Declining infected")
dfd2=merge(df2,dfdI,all=T)
ggplot(dfd2,aes(x=time,y=N))+geom_point(aes(colour=species))+facet_wrap(~status)


## leave em and mm constant and vary me and ee
alphaC.em=alpha.em[10]
alphaC.mm=alpha.mm[2]

# simulation 1
alpha.ee1=alpha.ee[2]
alpha.me1=alpha.me[2]

# Ev by itself
Ni.ev1=rep(NA,simtime)
Ni.ev1[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev1[t]=Ni.ev1[t-1]*germ.ev*lambda.ev/(1+alpha.ee1*germ.ev*Ni.ev1[t-1])
}

# Infected Mv population
NdI.mv1=rep(NA,simtime)
NdI.mv1[1]=N0.mv

# Healthy Ev population
NdH.ev1=rep(NA,simtime)
NdH.ev1[1]=Ni.ev1[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev1[t]=NdH.ev1[t-1]*germ.ev*lambda.ev/(1+alpha.ee1*germ.ev*NdH.ev1[t-1]+alphaC.em*germ.mv*NdI.mv1[t-1])
	NdI.mv1[t]=NdI.mv1[t-1]*germ.mv*lambdadI.mv[t]/(1+alphaC.mm*germ.mv*NdI.mv1[t-1]+alpha.me1*germ.ev*NdH.ev1[t-1])
	NdH.ev1[t]=ifelse(NdH.ev1[t]<1,0,NdH.ev1[t])
	NdI.mv1[t]=ifelse(NdI.mv1[t]<1,0,NdI.mv1[t])
}

# plot
dfdI1=data.frame(N=c(NdH.ev1,NdI.mv1),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),coefficients="low, low")
ggplot(dfdI1,aes(x=time,y=N))+geom_point(aes(colour=species))


# simulation 2
alpha.ee2=alpha.ee[2]
alpha.me2=alpha.me[8]

# Ev by itself
Ni.ev2=rep(NA,simtime)
Ni.ev2[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev2[t]=Ni.ev2[t-1]*germ.ev*lambda.ev/(1+alpha.ee2*germ.ev*Ni.ev2[t-1])
}

# Infected Mv population
NdI.mv2=rep(NA,simtime)
NdI.mv2[1]=N0.mv

# Healthy Ev population
NdH.ev2=rep(NA,simtime)
NdH.ev2[1]=Ni.ev2[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev2[t]=NdH.ev2[t-1]*germ.ev*lambda.ev/(1+alpha.ee2*germ.ev*NdH.ev2[t-1]+alphaC.em*germ.mv*NdI.mv2[t-1])
	NdI.mv2[t]=NdI.mv2[t-1]*germ.mv*lambdadI.mv[t]/(1+alphaC.mm*germ.mv*NdI.mv2[t-1]+alpha.me2*germ.ev*NdH.ev2[t-1])
	NdH.ev2[t]=ifelse(NdH.ev2[t]<1,0,NdH.ev2[t])
	NdI.mv2[t]=ifelse(NdI.mv2[t]<1,0,NdI.mv2[t])
}

# plot
dfdI2=data.frame(N=c(NdH.ev2,NdI.mv2),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),coefficients="low, high")
ggplot(dfdI2,aes(x=time,y=N))+geom_point(aes(colour=species))

# simulation 3
alpha.ee3=alpha.ee[8]
alpha.me3=alpha.me[8]

# Ev by itself
Ni.ev3=rep(NA,simtime)
Ni.ev3[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev3[t]=Ni.ev3[t-1]*germ.ev*lambda.ev/(1+alpha.ee3*germ.ev*Ni.ev3[t-1])
}

# Infected Mv population
NdI.mv3=rep(NA,simtime)
NdI.mv3[1]=N0.mv

# Healthy Ev population
NdH.ev3=rep(NA,simtime)
NdH.ev3[1]=Ni.ev3[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev3[t]=NdH.ev3[t-1]*germ.ev*lambda.ev/(1+alpha.ee3*germ.ev*NdH.ev3[t-1]+alphaC.em*germ.mv*NdI.mv3[t-1])
	NdI.mv3[t]=NdI.mv3[t-1]*germ.mv*lambdadI.mv[t]/(1+alphaC.mm*germ.mv*NdI.mv3[t-1]+alpha.me3*germ.ev*NdH.ev3[t-1])
	NdH.ev3[t]=ifelse(NdH.ev3[t]<1,0,NdH.ev3[t])
	NdI.mv3[t]=ifelse(NdI.mv3[t]<1,0,NdI.mv3[t])
}

# plot
dfdI3=data.frame(N=c(NdH.ev3,NdI.mv3),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),coefficients="high, high")
ggplot(dfdI3,aes(x=time,y=N))+geom_point(aes(colour=species))

# simulation 4
alpha.ee4=alpha.ee[5]
alpha.me4=alpha.me[2]

# Ev by itself
Ni.ev4=rep(NA,simtime)
Ni.ev4[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev4[t]=Ni.ev4[t-1]*germ.ev*lambda.ev/(1+alpha.ee4*germ.ev*Ni.ev4[t-1])
}

# Infected Mv population
NdI.mv4=rep(NA,simtime)
NdI.mv4[1]=N0.mv

# Healthy Ev population
NdH.ev4=rep(NA,simtime)
NdH.ev4[1]=Ni.ev4[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev4[t]=NdH.ev4[t-1]*germ.ev*lambda.ev/(1+alpha.ee4*germ.ev*NdH.ev4[t-1]+alphaC.em*germ.mv*NdI.mv4[t-1])
	NdI.mv4[t]=NdI.mv4[t-1]*germ.mv*lambdadI.mv[t]/(1+alphaC.mm*germ.mv*NdI.mv4[t-1]+alpha.me4*germ.ev*NdH.ev4[t-1])
	NdH.ev4[t]=ifelse(NdH.ev4[t]<1,0,NdH.ev4[t])
	NdI.mv4[t]=ifelse(NdI.mv4[t]<1,0,NdI.mv4[t])
}

# plot
dfdI4=data.frame(N=c(NdH.ev4,NdI.mv4),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),intraC=5)
ggplot(dfdI4,aes(x=time,y=N))+geom_point(aes(colour=species))

# simulation 5
alpha.ee5=alpha.ee[6]
alpha.me5=alpha.me[2]

# Ev by itself
Ni.ev5=rep(NA,simtime)
Ni.ev5[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev5[t]=Ni.ev5[t-1]*germ.ev*lambda.ev/(1+alpha.ee5*germ.ev*Ni.ev5[t-1])
}

# Infected Mv population
NdI.mv5=rep(NA,simtime)
NdI.mv5[1]=N0.mv

# Healthy Ev population
NdH.ev5=rep(NA,simtime)
NdH.ev5[1]=Ni.ev5[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev5[t]=NdH.ev5[t-1]*germ.ev*lambda.ev/(1+alpha.ee5*germ.ev*NdH.ev5[t-1]+alphaC.em*germ.mv*NdI.mv5[t-1])
	NdI.mv5[t]=NdI.mv5[t-1]*germ.mv*lambdadI.mv[t]/(1+alphaC.mm*germ.mv*NdI.mv5[t-1]+alpha.me5*germ.ev*NdH.ev5[t-1])
	NdH.ev5[t]=ifelse(NdH.ev5[t]<1,0,NdH.ev5[t])
	NdI.mv5[t]=ifelse(NdI.mv5[t]<1,0,NdI.mv5[t])
}

# plot
dfdI5=data.frame(N=c(NdH.ev5,NdI.mv5),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),intraC=6)
ggplot(dfdI5,aes(x=time,y=N))+geom_point(aes(colour=species))

# simulation 6
alpha.ee6=alpha.ee[7]
alpha.me6=alpha.me[2]

# Ev by itself
Ni.ev6=rep(NA,simtime)
Ni.ev6[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev6[t]=Ni.ev6[t-1]*germ.ev*lambda.ev/(1+alpha.ee6*germ.ev*Ni.ev6[t-1])
}

# Infected Mv population
NdI.mv6=rep(NA,simtime)
NdI.mv6[1]=N0.mv

# Healthy Ev population
NdH.ev6=rep(NA,simtime)
NdH.ev6[1]=Ni.ev6[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev6[t]=NdH.ev6[t-1]*germ.ev*lambda.ev/(1+alpha.ee6*germ.ev*NdH.ev6[t-1]+alphaC.em*germ.mv*NdI.mv6[t-1])
	NdI.mv6[t]=NdI.mv6[t-1]*germ.mv*lambdadI.mv[t]/(1+alphaC.mm*germ.mv*NdI.mv6[t-1]+alpha.me6*germ.ev*NdH.ev6[t-1])
	NdH.ev6[t]=ifelse(NdH.ev6[t]<1,0,NdH.ev6[t])
	NdI.mv6[t]=ifelse(NdI.mv6[t]<1,0,NdI.mv6[t])
}

# plot
dfdI6=data.frame(N=c(NdH.ev6,NdI.mv6),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),intraC=7)
ggplot(dfdI6,aes(x=time,y=N))+geom_point(aes(colour=species))


# simulation 7
alpha.ee7=alpha.ee[8]
alpha.me7=alpha.me[2]

# Ev by itself
Ni.ev7=rep(NA,simtime)
Ni.ev7[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev7[t]=Ni.ev7[t-1]*germ.ev*lambda.ev/(1+alpha.ee7*germ.ev*Ni.ev7[t-1])
}

# Infected Mv population
NdI.mv7=rep(NA,simtime)
NdI.mv7[1]=N0.mv

# Healthy Ev population
NdH.ev7=rep(NA,simtime)
NdH.ev7[1]=Ni.ev7[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev7[t]=NdH.ev7[t-1]*germ.ev*lambda.ev/(1+alpha.ee7*germ.ev*NdH.ev7[t-1]+alphaC.em*germ.mv*NdI.mv7[t-1])
	NdI.mv7[t]=NdI.mv7[t-1]*germ.mv*lambdadI.mv[t]/(1+alphaC.mm*germ.mv*NdI.mv7[t-1]+alpha.me7*germ.ev*NdH.ev7[t-1])
	NdH.ev7[t]=ifelse(NdH.ev7[t]<1,0,NdH.ev7[t])
	NdI.mv7[t]=ifelse(NdI.mv7[t]<1,0,NdI.mv7[t])
}

# plot
dfdI7=data.frame(N=c(NdH.ev7,NdI.mv7),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),intraC=8)
ggplot(dfdI7,aes(x=time,y=N))+geom_point(aes(colour=species))

# simulation 4b, no infection

# Ev by itself
Ni.ev4b=rep(NA,simtime)
Ni.ev4b[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	Ni.ev4b[t]=Ni.ev4b[t-1]*germ.ev*lambda.ev/(1+alpha.ee4*germ.ev*Ni.ev4b[t-1])
}

# Infected Mv population
NdH.mv4b=rep(NA,simtime)
NdH.mv4b[1]=N0.mv

# Healthy Ev population
NdH.ev4b=rep(NA,simtime)
NdH.ev4b[1]=Ni.ev4b[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev4b[t]=NdH.ev4b[t-1]*germ.ev*lambda.ev/(1+alpha.ee4*germ.ev*NdH.ev4b[t-1]+alphaC.em*germ.mv*NdH.mv4b[t-1])
	NdH.mv4b[t]=NdH.mv4b[t-1]*germ.mv*lambda.mv/(1+alphaC.mm*germ.mv*NdH.mv4b[t-1]+alpha.me4*germ.ev*NdH.ev4b[t-1])
	NdH.ev4b[t]=ifelse(NdH.ev4b[t]<1,0,NdH.ev4b[t])
	NdH.mv4b[t]=ifelse(NdH.mv4b[t]<1,0,NdH.mv4b[t])
}

# plot
dfdI4b=data.frame(N=c(NdH.ev4b,NdH.mv4b),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),intraC=3)
ggplot(dfdI4b,aes(x=time,y=N))+geom_point(aes(colour=species))


# merge 4 through 7
dfdIEE=rbind(dfdI4b,dfdI4,dfdI5,dfdI6,dfdI7)

# plot
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("MvEvPAIDModel_103018.pdf",width=9,height=5.5)
ggplot(dfdIEE,aes(x=time,y=N))+geom_line(aes(colour=species),size=1.5)+facet_wrap(~as.factor(intraC),nrow=1)+scale_colour_manual(values=c("tomato1","cyan4"),name="Species")+theme_bw()+theme(legend.position=c(0.9,0.8),strip.text=element_blank(),axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+xlab("Days")+ylab("Population size")
dev.off()