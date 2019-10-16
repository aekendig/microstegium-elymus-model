# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(data.table)
library(plotly)
library(tidyr)
library(cowplot)
library(popbio)

# plotting parameters
axisText=24
axisTitle=26
legendText=20
legendTitle=22

## Parameters
# see ModelBuilding powerpoint for details
# subscripts follow "."
# p = perennial adult (at least 1 year old)
# s = perennial seedling (germinated that spring)
# a = annual
# L = annual litter

# perennial adult survival
m.p=0.95

# dormant seed survival
s.a=1-0.26
s.s = 1

# germination in the absence of competition
gamma.a=0.86
gamma.s=0.8

# reduction in germination due to litter
alpha.aL=0.13
alpha.sL=0.027

# proportion of Mv biomass that becomes litter the next year
b=0.56

# survival to reproduction
h.p=0.95
h.a=0.90
h.s.mean=0.46
h.s.int=-0.8258
h.s.slope=0.2859

# seed production rates
lambda.a=6500 
lambda.p=435

# seed production competition parameters
alpha.aa=mean(c(0.015,0.001))
alpha.as=0.0037
alpha.sa=0.091
alpha.ps=0
alpha.pp=0.06

# convert seed production to biomass
c.a=30/6500

# decrease seed production due to infection
tolMax=0.31
tol.a=0.6
tol.p=0.95

# years
simtime.p=100
simtime.i=300

# initial numbers
N0.a=1
N0.s=0
N0.p=1

# parameter list
params=c(m.p,s.a,s.s,gamma.a,gamma.s,alpha.aL,alpha.sL,b,h.p,h.a,lambda.a,lambda.p,alpha.aa,alpha.as,alpha.sa,alpha.pp,c.a,tol.a,tol.p)


## Population Model

simFun=function(params,h.s.fac,invader,N0.a,N0.s,N0.p,L0){
	
	# define parameters
	m.p=params[1]
	s.a=params[2]
	s.s=params[3]
	gamma.a=params[4]
	gamma.s=params[5]
	alpha.aL=params[6]
	alpha.sL=params[7]
	b=params[8]
	h.p=params[9]
	h.a=params[10]
	lambda.a=params[11]
	lambda.p=params[12]
	alpha.aa=params[13]
	alpha.as=params[14]
	alpha.sa=params[15]
	alpha.pp=params[16]
	c.a=params[17]
	tol.a=params[18]
	tol.p=params[19]
	
	# calculate parameters
	lambda.s=lambda.p/10
	alpha.ap=alpha.as*10
	alpha.ss=alpha.as
	alpha.sp=alpha.ap
	alpha.pa=alpha.sa/10
	lam.a=tol.a*lambda.a
	lam.p=tol.p*lambda.p
	lam.s=tol.p*lambda.s	
	
	# functions
  h.s_fun=function(x){
    h.out=exp(h.s.int+h.s.slope*x)/(1+exp(h.s.int+h.s.slope*x))
    return(h.out)
  }

	# initialize populations
	N.a=rep(NA,simtime.i)
	N.s=rep(NA,simtime.i)
	N.p=rep(NA,simtime.i)
	L=rep(NA,simtime.i)

	N.a[1]=N0.a
	N.s[1]=N0.s
	N.p[1]=N0.p
	L[1]=L0

	# simulate population dynamics
	for(t in 1:(simtime.i-1)){	
		
		# composite paramters
	  	h.s=ifelse(h.s.fac==T,h.s_fun(N.s[t]),h.s.mean)
	  	h.s=ifelse(is.na(h.s)&h.s.fac==T,1,h.s)
	  	g.s=gamma.s/(1+alpha.sL*L[t])
	  	g.a=gamma.a/(1+alpha.aL*L[t])
	  	f.s=lam.s/(1+alpha.ss*g.s*h.s*N.s[t]+alpha.sp*m.p*N.p[t]+alpha.sa*g.a*h.a*N.a[t])
	  	f.p=lam.p/(1+alpha.ps*g.s*h.s*N.s[t]+alpha.pp*m.p*N.p[t]+alpha.pa*g.a*h.a*N.a[t])
	  	f.a=lam.a/(1+alpha.as*g.s*h.s*N.s[t]+alpha.ap*m.p*N.p[t]+alpha.aa*g.a*h.a*N.a[t])

		# population size
		N.s[t+1]=s.s*(1-g.s)*N.s[t]+g.s*h.s*f.s*N.s[t]+m.p*f.p*N.p[t]
		N.p[t+1]=m.p*N.p[t]+g.s*h.s*N.s[t]
		N.a[t+1]=s.a*(1-g.a)*N.a[t]+g.a*h.a*f.a*N.a[t]
		L[t+1]=b*c.a*g.a*h.a*f.a*N.a[t]+b*L[t]
		
		# correct to prevent negative numbers
		N.s[t+1]=ifelse(N.s[t+1]<1,0,N.s[t+1])
		N.p[t+1]=ifelse(N.p[t+1]<1,0,N.p[t+1])
		N.a[t+1]=ifelse(N.a[t+1]<1,0,N.a[t+1])
		L[t+1]=ifelse(L[t+1]<0,0,L[t+1])
	}

	# perennial transition matrix
	pmat=matrix(c(s.s*(1-g.s)+g.s*h.s*f.s,m.p*f.p,g.s*h.s,m.p),nrow=2,byrow=T)
	grwr.p=eigen.analysis(pmat)$lambda1
	
	# annual transition matrix
	amat=matrix(c(s.a*(1-g.a)+g.a*h.a*f.a,0,b*c.a*g.a*h.a*f.a,b),nrow=2,byrow=T)
	grwr.a=eigen.analysis(amat)$lambda1
	
	# export grwr
	grwr=ifelse(invader=="P",grwr.p,grwr.a)
	
	# save data
	dfN=data.frame(time=rep(1:simtime.i,4),N=c(N.s,N.p,N.a,L),species=rep(c("Elymus seedling","Elymus adult","Microstegium","Microstegium litter"),each=simtime.i))
	dfN$grwr=grwr
	
	# return
	return(dfN)
}


## Run simulation with default values

# Perennial invasion
pinv=simFun(params,h.s.fac=F,invader="P",N0.a=1,N0.s=0,N0.p=0,L0=0)

# Perennial grwr
unique(pinv$grwr) #1.2
log(unique(pinv$grwr)) #0.18

# Check simulation
ggplot(pinv)+geom_line(aes(x=time,y=N,colour=species),size=1)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")

tail(pinv)

# Annual invasion
ainv=simFun(params,h.s.fac=F,invader="A",N0.a=0,N0.s=0,N0.p=1,L0=0)

# Annual grwr
unique(ainv$grwr) #1.4
log(unique(ainv$grwr)) #0.36

# Check simulation
ggplot(ainv)+geom_line(aes(x=time,y=N,colour=species),size=1)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")

tail(ainv)