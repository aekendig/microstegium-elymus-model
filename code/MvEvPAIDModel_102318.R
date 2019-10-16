# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(geepack)
library(doBy)
library(lme4)
library(Hmisc)
library(viridis)

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
(lambda.mv=2500/10*12.5/20)
lambda.ev=435

# initial numbers
N0.mv=50
N0.ev=50

# intraspecific competition
alpha.mm=seq(0,1,by=0.1)
alpha.ee=seq(0,1,by=0.1)

# interspecific competition
alpha.me=seq(0,1,by=0.1)
alpha.em=seq(0,1,by=0.1)

# years
simtime=100

# poplations
N.mv=rep(NA,simtime)
N.ev=rep(NA,simtime)

N.mv[1]=N0.mv
N.ev[1]=N0.ev

# intialize Ev population
for(t in 2:simtime){
	N.ev[t]=N.ev[t-1]*germ.ev*lambda.ev/(1+alpha.ee[8]*germ.ev*N.ev[t-1])
}

# plot
df.ev=data.frame(N=N.ev,time=1:simtime)
ggplot(df.ev,aes(x=time,y=N))+geom_point()