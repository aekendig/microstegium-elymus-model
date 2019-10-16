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
(germ.a=mean(mvDat$germ.mv))
(germ.s=mean(evDat$germ.ev))

# seed production rates
(lambda.a=2500/10*12.5/20) # Wilson et al. 2015
lambda.p=435 # Stevens 1957

# intraspecific/group competition
alpha.aa=0.2 # annuals
alpha.ss=0.8 # perennial seedlings
alpha.pp=0.2 # perennial adults

# interspecific/group competition
alpha.as=0.4 # p seedlings on annuals
alpha.ap=0.6 # p adults on annuals
alpha.sa=2 # annuals on p seedlings
alpha.sp=1 # p adults on p seedlings
alpha.pa=0.2 # annuals on p adults
alpha.ps=0.2 # p seedlings on p adults

# adult perennial survival
s.p=0.94

# years
simtime=100

# initial numbers
N0.a=1
N0.s=1
N0.p=0


## Perennial starts off by itself

# initialize populations
Ni.s=rep(NA,simtime)
Ni.s[1]=N0.s
Ni.p=rep(NA,simtime)
Ni.p[1]=N0.p

# perennial population dynamics
for(t in 1:(simtime-1)){
	Ni.s[t+1]=Ni.p[t]*lambda.p/(1+alpha.ps*germ.s*Ni.s[t]+alpha.pp*Ni.p[t])
	Ni.p[t+1]=Ni.s[t]*germ.s/(1+alpha.ss*germ.s*Ni.s[t]+alpha.sp*Ni.p[t])+Ni.p[t]*s.p
}

# plot
dfNi=data.frame(time=1:simtime,N.s=Ni.s,N.p=Ni.p)
ggplot(dfNi)+geom_line(aes(x=time,y=N.s),colour="green")+geom_line(aes(x=time,y=N.p),colour="blue")


## Annual invades perennial

# initialize populations
Nh.a=rep(NA,simtime)
Nh.s=rep(NA,simtime)
Nh.p=rep(NA,simtime)

Nh.a[1]=N0.a
Nh.s[1]=Ni.s[simtime]
Nh.p[1]=Ni.p[simtime]

# introduce annual population
for(t in 1:(simtime-1)){
	Nh.s[t+1]=Nh.p[t]*lambda.p/(1+alpha.pa*germ.a*Nh.a[t]+alpha.ps*germ.s*Nh.s[t]+alpha.pp*Nh.p[t])
	Nh.p[t+1]=Nh.s[t]*germ.s/(1+alpha.sa*germ.a*Nh.a[t]+alpha.ss*germ.s*Nh.s[t]+alpha.sp*Nh.p[t])+Nh.p[t]*s.p
	Nh.a[t+1]=Nh.a[t]*germ.a*lambda.a/(1+alpha.aa*germ.a*Nh.a[t]+alpha.as*germ.s*Nh.s[t]+alpha.ap*Nh.p[t])
}

# plot
dfN=data.frame(time=1:simtime,N.s=Nh.s,N.p=Nh.p,N.a=Nh.a)
ggplot(dfN)+geom_line(aes(x=time,y=N.s),colour="green")+geom_line(aes(x=time,y=N.p),colour="blue")+geom_line(aes(x=time,y=N.a),colour="red")


## Add infection to system

# initialize populations
N.a=rep(NA,simtime)
N.s=rep(NA,simtime)
N.p=rep(NA,simtime)

N.a[1]=N0.a
N.s[1]=Ni.s[simtime]
N.p[1]=Ni.p[simtime]

# SIR model
SIRModel=function(
g.A,g.P,g.Q,
delta.AA,delta.AP,delta.AQ,
delta.PA,delta.PP,delta.PQ,
delta.QA,delta.QP,delta.QQ,
beta.AA,beta.AP,beta.AQ,
beta.PA,beta.PP,beta.PQ,
beta.QA,beta.QP,beta.QQ,
A0.s,P0.s,Q0.s,
A0.i,P0.i,Q0.i,
mu.As,mu.Ai,mu.Ps,
mu.Pi,mu.Qs,mu.Qi){
	
	parms=list(BP=BP, qP=qP, vP=vP, BR=BR, qR=qR, vR=vR, vC=vC, Pinit=Pinit, Rinit=Rinit, Cinit=Cinit, simtime=simtime)
	
	mymodel=with(as.list(parms),function(t,x,parms){
		
		A.s=x["A.s"]
		A.i=x["A.i"]
		A=A.s+A.i
		P.s=x["P.s"]
		P.i=x["P.i"]
		P=P.s+P.i
		Q.s=x["Q.s"]
		Q.i=x["Q.i"]
		Q=Q.s+Q.i

### START HERE ###
		A.s.dot = g.A*A.s/(1+delta.AA*A+delta.AP*P+delta.AQ*Q)-beta.AA*A.s*A.i-
		Rdot = BR*(R+qR*C)*S - BP*(exp(logP)+qP*C)*R - vR*R
		Cdot = BP*(exp(logP)+qP*C)*R + BR*(R+qR*C)*exp(logP) - vC*C
		
		list(c(logPdot, Rdot, Cdot))
	})
	
	xstart=c(logP=log(Pinit), R=Rinit, C=Cinit)
	
	times=seq(0,simtime,length=simtime)
	
	out=as.data.frame(lsoda(xstart, times, mymodel, parms, hmax=20))
	out$S=with(out,1-exp(logP)-R-C)
	
	return(out)
}

# introduce annual population and infection
for(t in 1:(simtime-1)){
	
	
	
	
	
	N.s[t+1]=N.p[t]*lambda.p/(1+alpha.pa*germ.a*N.a[t]+alpha.ps*germ.s*N.s[t]+alpha.pp*N.p[t])
	N.p[t+1]=N.s[t]*germ.s/(1+alpha.sa*germ.a*N.a[t]+alpha.ss*germ.s*N.s[t]+alpha.sp*N.p[t])+N.p[t]*s.p
	N.a[t+1]=N.a[t]*germ.a*lambda.a/(1+alpha.aa*germ.a*N.a[t]+alpha.as*germ.s*N.s[t]+alpha.ap*N.p[t])
}

# plot
dfN=data.frame(time=1:simtime,N.s=N.s,N.p=N.p,N.a=N.a)
ggplot(dfN)+geom_line(aes(x=time,y=N.s),colour="green")+geom_line(aes(x=time,y=N.p),colour="blue")+geom_line(aes(x=time,y=N.a),colour="red")
