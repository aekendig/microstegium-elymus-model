# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(data.table)
library(plotly)
library(tidyr)
library(cowplot)

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
virMax=0.31
vir.a=0.6
vir.p=0.95

# years
simtime.p=100
simtime.i=300

# initial numbers
N0.a=1
N0.s=0
N0.p=1

# parameter list
params=c(m.p,s.a,s.s,gamma.a,gamma.s,alpha.aL,alpha.sL,b,h.p,h.a,lambda.a,lambda.p,alpha.aa,alpha.as,alpha.sa,alpha.pp,c.a,vir.a,vir.p)


## Population Model

simFun=function(params,h.s.comp,invader,N0.a,N0.s,N0.p,L0){
	
	# calculate parameters
	lambda.s=lambda.p/10
	alpha.ap=alpha.as*10
	alpha.ss=alpha.as
	alpha.sp=alpha.ap
	alpha.pa=alpha.sa/10
	lam.a=vir.a*lambda.a
	lam.p=vir.p*lambda.p
	lam.s=vir.p*lambda.s	
	
	# functions
  h.s_fun=function(x){
    h.out=exp(h.s.int+h.s.slope*x)/(1+exp(h.s.int+h.s.slope*x))
    return(h.out)
  }

	# initialize populations
	Nh.a=rep(NA,simtime.i)
	Nh.s=rep(NA,simtime.i)
	Nh.p=rep(NA,simtime.i)
	Lh=rep(NA,simtime.i)

	Nh.a[1]=N0.a
	Nh.s[1]=N0.s
	Nh.p[1]=N0.p
	Lh[1]=L0

	# simulate population dynamics
	for(t in 1:(simtime.i-1)){	
		
		# composite paramters
	  	h.s=ifelse(h.s.comp==T,h.s_fun(Nh.s[t]),h.s.mean)
	  	h.s=ifelse(is.na(h.s)&h.s.comp==T,1,h.s)
	  	g.s=gamma.s/(1+alpha.sL*Lh[t])
	  	g.a=gamma.a/(1+alpha.aL*Lh[t])
	  	f.s=lam.s/(1+alpha.ss*g.s*h.s*Nh.s[t]+alpha.sp*m.p*Nh.p[t]+alpha.sa*g.a*h.a*Nh.a[t])
	  	f.p=lam.p/(1+alpha.ps*g.s*h.s*Nh.s[t]+alpha.pp*m.p*Nh.p[t]+alpha.pa*g.a*h.a*Nh.a[t])
	  	f.a=lam.a/(1+alpha.as*g.s*h.s*Nh.s[t]+alpha.ap*m.p*Nh.p[t]+alpha.aa*g.a*h.a*Nh.a[t])

		# population size
		Nh.s[t+1]=s.s*(1-g.s)*Nh.s[t]+g.s*h.s*f.s*Nh.s[t]+m.p*f.p*Nh.p[t]
		Nh.p[t+1]=m.p*Nh.p[t]+g.s*h.s*Nh.s[t]
		Nh.a[t+1]=s.a*(1-g.a)*Nh.a[t]+g.a*h.a*f.a*Nh.a[t]
		Lh[t+1]=b*c.a*g.a*h.a*f.a*Nh.a[t]+b*Lh[t]
		
		# correct to prevent negative numbers
		Nh.s[t+1]=ifelse(Nh.s[t+1]<1,0,Nh.s[t+1])
		Nh.p[t+1]=ifelse(Nh.p[t+1]<1,0,Nh.p[t+1])
		Nh.a[t+1]=ifelse(Nh.a[t+1]<1,0,Nh.a[t+1])
		Bh.a[t+1]=ifelse(Bh.a[t+1]<0,0,Bh.a[t+1])
	}

	# perennial transition matrix
	pmat=matrix(c(s.s*(1-g.s)+g.s*h.s*f.s,m.p*f.p,g.s*h.s,m.p),nrow=2,byrow=T)
	grwr.p=eigen(pmat)
	
	# annual transition matrix
	amat=matrix(c(s.a*(1-g.a)+g.a*h.a*f.a,0,b*c.a*g.a*h.a*f.a,b),nrow=2,byrow=T)
	grwr.a=eigen(amat)
	
	# export grwr
	grwr=ifelse(invader=="P",grwr.p,grwr.a)
	
	
	# save data
	dfNh=data.frame(time=rep(1:simtime.i,4),N=c(Nh.s,Nh.p,Nh.a,Bh.a),species=rep(c("Elymus seedling","Elymus adult","Microstegium","Microstegium biomass"),each=simtime.i))
	dfNh$scenario="Invasion"
	
	# combine scenarios
	dfTot=rbind(dfNi,dfNh)
	dfTot$aa=alpha.aa
	dfTot$pp=alpha.pp
	dfTot$vir.p=vir.p
	
	# return
	return(dfTot)
}

### TESTING SIMULATION IN GRWR SCRIPT, BELOW IS FROM PREVIOUS VERSION ###


## Run simulation with default values

outdf=compFun(params=params,h.s.comp=F)
outdf$scenario=factor(outdf$scenario,levels=c("Pre-invasion","Invasion"))

ggplot(outdf)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_wrap(~scenario,scales="free_x")+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")


### STILL WRITING BELOW FUNCTION ###

## Sensitivity analysis function
sensFunc=function(npar,P,pp,dp,p,inits,path_init,model,parm_names,min.y,max.y){

	# Calculate sensitivities for parameters and initial values
	mabs=rep(NA,npar) 
	msqr=rep(NA,npar) 
	ms=rep(NA,npar)
	sens=rep(NA,npar)


	for (i in 1:npar){
		# Add small perturbation to a parameter
		pstar=paramsUp[i]  
		# Replace parameter with perturbed value     
  		p[i]=pstar
  		# Solve ODE with perturbed value
  		outstar=compFun(params=p,h.s.comp=F) 
  		# Virulence with perturbation
  		virulence_star = (outstar_healthy[,4] - outstar[,4])/outstar_healthy[,4]
  		Pstar=virulence_star[eval_interval]
  		Rel = pp[i]/P
  		# Calculate sensitivity
  		sens[,i]=(Pstar-P)/(pstar-pp[i])*pp[i]/P  
  		# Replace perturbed value with original
  		p[i]=pp[i]              
  		linetype=ifelse(i<9,1,ifelse(i>16,3,2))
  		# Plot sensitivity
  		lines(times[eval_interval],sens[,i],col=i,lty=linetype,lwd=2) 
  		
  		# Mean value of sensitivity
  		ms[i]=mean(sens[,i]) 
  		# Mean absolute value of sensitivity
  		mabs[i]=mean(abs(sens[,i])) 
  		# Mean square of sensitivity
  		msqr[i]=sqrt(sum(sens[,i]^2)/length(sens[,i])) 
	}

	pnames=parm_names
	legend("topleft",pnames,col=c(1:length(pnames)),cex=1,lwd=2,lty=c(rep(1,8),rep(2,8),rep(3,8)),ncol=6)
	
	# Mean sensitivity plots
		barplot(ms, main = "Mean sensitivity", xlab = "parameter", ylab = "virulence (mass lost)", names.arg = pnames,cex.names=0.5)
	barplot(mabs, main = "Mean absolute value\n of sensitivity", xlab = "parameter", ylab = "virulence (mass lost)", names.arg = pnames,cex.names=0.5)
	barplot(msqr, main = "Mean square of sensitivity", xlab = "parameter", ylab = "virulence (mass lost)", names.arg = pnames,cex.names=0.5)
	
	# Output dataframe
	dat.out=data.frame(parms=pnames,ms=ms,mabs=mabs,msqr=msqr)
}



## Adjust parameter values

paramsUp=params