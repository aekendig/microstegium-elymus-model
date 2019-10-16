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

# parameter names
paramNames=c("m.p","s.a","s.s","gamma.a","gamma.s","alpha.aL","alpha.sL","b","h.p","h.a","lambda.a","lambda.p","alpha.aa","alpha.as","alpha.sa","alpha.pp","c.a","tol.a","tol.p")


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


## Sensitivity analysis function
sensFun=function(sensVal,paramNames,params,invader,N0.a,N0.p){

	# run model with default parameters
out_def=simFun(params=params,h.s.fac=F,invader=invader,N0.a=N0.a,N0.s=0,N0.p=N0.p,L0=0)

	# save default output
	grwr_def=log(unique(out_def$grwr))
	
	# increase parameters by 5%
	paramsAdj=params*sensVal
	
	# number of parameters
	npar=length(params)
	
	# output vectors
	sens=rep(NA,npar)
	grwr_star=rep(NA,npar)

	# cycle through values and adjust each
	for (i in 1:npar){
		
		# create temporary parameter string
		p=params
		
		# add small perturbation to a parameter
		p[i]=paramsAdj[i]  
  		
  		# run model with perturbed value
  		out_star=simFun(params=p,h.s.fac=F,invader=invader,N0.a=N0.a,N0.s=0,N0.p=N0.p,L0=0) 
  		
  		# grwr with perturbed values
  		grwr_star[i]=log(unique(out_star$grwr))
  		
  		# calculate sensitivity
  		sens[i]=(grwr_star[i]-grwr_def)/(p[i]-params[i])*params[i]/grwr_def  
	}

	# output data frame
	sensdf=data.frame(param=paramNames,value=params,valueAdj=paramsAdj,grwr=rep(grwr_def,npar),grwrAdj=grwr_star,sensitivity=sens)
	return(sensdf)
}


## Sensitivity analysis

# Increase by 5%

# Perennial invasion
pinv=sensFun(sensVal=1.05,paramNames=paramNames,params=params,invader="P",N0.a=1,N0.p=0)
pinv

# Annual invasion
ainv=sensFun(sensVal=1.05,paramNames=paramNames,params=params,invader="A",N0.a=0,N0.p=1)
ainv

# Decrease by 5%

# Perennial invasion
pinv2=sensFun(sensVal=0.95,paramNames=paramNames,params=params,invader="P",N0.a=1,N0.p=0)
pinv2

# Annual invasion
ainv2=sensFun(sensVal=0.95,paramNames=paramNames,params=params,invader="A",N0.a=0,N0.p=1)
ainv2

# Combine dataframes
pinv$invader="perennial log(GRWR)"
ainv$invader="annual log(GRWR)"
pinv2$invader="perennial log(GRWR)"
ainv2$invader="annual log(GRWR)"

pinv$adjustment="increase 5%"
ainv$adjustment="increase 5%"
pinv2$adjustment="decrease 5%"
ainv2$adjustment="decrease 5%"

invdf=rbind(pinv,ainv,pinv2,ainv2)

# Absolute value of sensitivity
invdf$sensitivityAbs=abs(invdf$sensitivity)

# Plot
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("FullSensitivityAnalysis_PlantDiseaseModel_021219.pdf")
ggplot(invdf,aes(x=param,y=sensitivity))+geom_bar(stat="identity")+facet_grid(adjustment~invader)+geom_hline(yintercept=0,size=0.3)+theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("SensitivityAnalysis_PlantDiseaseModel_021219.pdf")
ggplot(filter(invdf,adjustment=="increase 5%"),aes(x=param,y=sensitivity))+geom_bar(stat="identity")+facet_wrap(~invader)+geom_hline(yintercept=0,size=0.3)+theme(axis.text.x=element_text(angle=45,hjust=1,size=10),strip.text=element_text(size=10))+xlab("parameter")
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
write.csv(invdf,"SensitivityAnalysis_PlantDiseaseModel_021219.csv",row.names=F)


## Repeat with no infection

paramsH=c(m.p,s.a,s.s,gamma.a,gamma.s,alpha.aL,alpha.sL,b,h.p,h.a,lambda.a,lambda.p,alpha.aa,alpha.as,alpha.sa,alpha.pp,c.a,1,1)

# Increase by 5%

# Perennial invasion
pinvH=sensFun(sensVal=1.05,paramNames=paramNames,params=paramsH,invader="P",N0.a=1,N0.p=0)
pinvH

# Annual invasion
ainvH=sensFun(sensVal=1.05,paramNames=paramNames,params=paramsH,invader="A",N0.a=0,N0.p=1)
ainvH

# Combine dataframes
pinvH$invader="perennial log(GRWR)"
ainvH$invader="annual log(GRWR)"

invHdf=rbind(pinvH,ainvH)

# Absolute value of sensitivity
invHdf$sensitivityAbs=abs(invHdf$sensitivity)

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("SensitivityAnalysis_HealthyPlantModel_021219.pdf")
ggplot(filter(invHdf,param!="tol.a"&param!="tol.p"),aes(x=param,y=sensitivity))+geom_bar(stat="identity")+facet_wrap(~invader)+geom_hline(yintercept=0,size=0.3)+theme(axis.text.x=element_text(angle=45,hjust=1,size=10),strip.text=element_text(size=10))+xlab("parameter")
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
write.csv(filter(invHdf,param!="tol.a"&param!="tol.p"),"SensitivityAnalysis_HealthyPlantModel_021219.csv",row.names=F)