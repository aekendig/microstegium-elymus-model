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
alpha.mm=seq(0,5,by=0.5)
alpha.ee=seq(0,5,by=0.5)

# interspecific competition
alpha.me=seq(0,5,by=0.5)
alpha.em=seq(0,5,by=0.5)

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
	N.ev[t]=N.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[8]*germ.ev*N.ev[t-1]+alpha.em[8]*germ.mv*N.mv[t-1])
	N.mv[t]=N.mv[t-1]*germ.mv*lambda.mv/(1+alpha.mm[2]*germ.mv*N.mv[t-1]+alpha.me[2]*germ.ev*N.ev[t-1])
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
	NH.ev[t]=NH.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[8]*germ.ev*NH.ev[t-1]+alpha.em[8]*germ.mv*NI.mv[t-1])
	NI.mv[t]=NI.mv[t-1]*germ.mv*lambdaI.mv/(1+alpha.mm[2]*germ.mv*NI.mv[t-1]+alpha.me[2]*germ.ev*NH.ev[t-1])
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
accTime=10
lambdadI.mv=c(seq(lambda.mv,lambdaI.mv,length=accTime),rep(lambdaI.mv,simtime-accTime))

# Infected Mv population
NdI.mv=rep(NA,simtime)
NdI.mv[1]=N0.mv

# Healthy Ev population
NdH.ev=rep(NA,simtime)
NdH.ev[1]=Ni.ev[simtime]

# introduce infected Mv population
for(t in 2:simtime){
	NdH.ev[t]=NdH.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[8]*germ.ev*NdH.ev[t-1]+alpha.em[8]*germ.mv*NdI.mv[t-1])
	NdI.mv[t]=NdI.mv[t-1]*germ.mv*lambdadI.mv[t]/(1+alpha.mm[2]*germ.mv*NdI.mv[t-1]+alpha.me[2]*germ.ev*NdH.ev[t-1])
}

# plot
dfdI=data.frame(N=c(NdH.ev,NdI.mv),time=rep(1:simtime,2),species=rep(c("Ev","Mv"),each=simtime),status="Declining infected")
dfd2=merge(df2,dfdI,all=T)
ggplot(dfd2,aes(x=time,y=N))+geom_point(aes(colour=species))+facet_wrap(~status)

