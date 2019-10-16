## need to update model to match model building powerpoint

# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(data.table)
library(plotly)

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
lambda.s=lambda.p/10

# seed production competition parameters
alpha.as=alpha.ss=alpha.ps=0
alpha.pa=0
#alpha.aa
#alpha.sa
#alpha.ap
#alpha.sp
#alpha.pp

# convert seed production to biomass
c.a=2.3e-4

# years
simtime=100

# initial numbers
N0.a=1
N0.s=0
N0.p=1

## Save N dynamics as a function and run across a range of intraspecific competition coefficients

# simulation function
compFun=function(alpha.aa,alpha.sa,alpha.ap,alpha.pp,h.s.comp){
	
	# functions
  h.s_fun=function(x){
    h.out=exp(h.s.int+h.s.slope*x)/(1+exp(h.s.int+h.s.slope*x))
  }
  
  # parameters
	alpha.sp=alpha.pp*0.1
	
	# perennial starts off by itself

	# initialize populations
	Ni.s=rep(NA,simtime)
	Ni.s[1]=N0.s
	Ni.p=rep(NA,simtime)
	Ni.p[1]=N0.p

	# perennial population dynamics
	for(t in 1:(simtime-1)){
	  
	  # composite parameters
	  h.s=ifelse(h.s.comp==T,h.s_fun(Ni.s[t]),h.s.mean)
	  g.s=gamma.s
	  f.s=lambda.s/(1+alpha.ss*g.s*h.s*Ni.s[t]+alpha.sp*m.p*Ni.p[t])
	  f.p=lambda.p/(1+alpha.ps*g.s*h.s*Ni.s[t]+alpha.pp*m.p*Ni.p[t])
	  
	  # population size
		Ni.s[t+1]=s.s*(1-g.s)*Ni.s[t]+g.s*h.s*f.s*Ni.s[t]+m.p*f.p*Ni.p[t]
		Ni.p[t+1]=m.p*Ni.p[t]+g.s*h.s*Ni.s[t]
		
		# correct to prevent negative numbers
		Ni.s[t+1]=ifelse(Ni.s[t+1]<1,0,Ni.s[t+1])
		Ni.p[t+1]=ifelse(Ni.p[t+1]<1,0,Ni.p[t+1])
	}

	# save data
	dfNi=data.frame(time=rep(1:simtime,2),N=c(Ni.s,Ni.p),species=rep(c("Elymus seedling","Elymus adult"),each=simtime))
	dfNi$scenario="Pre-invasion"

	# annual invades perennial

	# initialize populations
	Nh.a=rep(NA,simtime)
	Nh.s=rep(NA,simtime)
	Nh.p=rep(NA,simtime)
	Bh.a=rep(NA,simtime)

	Nh.a[1]=N0.a
	Nh.s[1]=Ni.s[simtime]
	Nh.p[1]=Ni.p[simtime]
	Bh.a[1]=0

	# introduce annual population
	for(t in 1:(simtime-1)){
	  
	  # composite paramters
	  h.s=ifelse(h.s.comp==T,h.s_fun(Nh.s[t]),h.s.mean)
	  Lh=b*Bh.a[t]
	  g.s=gamma.s/(1+alpha.sL*Lh)
	  g.a=gamma.a/(1+alpha.aL*Lh)
	  f.s=lambda.s/(1+alpha.ss*g.s*h.s*Nh.s[t]+alpha.sp*m.p*Nh.p[t]+alpha.sa*g.a*h.a*Nh.a[t])
	  f.p=lambda.p/(1+alpha.ps*g.s*h.s*Nh.s[t]+alpha.pp*m.p*Nh.p[t]+alpha.pa*g.a*h.a*Nh.a[t])
	  f.a=lambda.a/(1+alpha.as*g.s*h.s*Nh.s[t]+alpha.ap*m.p*Nh.p[t]+alpha.aa*g.a*h.a*Nh.a[t])
	  
	  # population size
		Nh.s[t+1]=s.s*(1-g.s)*Nh.s[t]+g.s*h.s*f.s*Nh.s[t]+m.p*f.p*Nh.p[t]
		Nh.p[t+1]=m.p*Nh.p[t]+g.s*h.s*Nh.s[t]
		Nh.a[t+1]=s.a*(1-g.a)*Nh.a[t]+g.a*h.a*f.a*Nh.a[t]
		Bh.a[t+1]=c.a*g.a*h.a*f.a*Nh.a[t]
		
		# correct to prevent negative numbers
		Nh.s[t+1]=ifelse(Nh.s[t+1]<1,0,Nh.s[t+1])
		Nh.p[t+1]=ifelse(Nh.p[t+1]<1,0,Nh.p[t+1])
		Nh.a[t+1]=ifelse(Nh.a[t+1]<1,0,Nh.a[t+1])
		Bh.a[t+1]=ifelse(Bh.a[t+1]<0,0,Bh.a[t+1])
	}

	# save data
	dfNh=data.frame(time=rep(1:simtime,4),N=c(Nh.s,Nh.p,Nh.a,Bh.a),species=rep(c("Elymus seedling","Elymus adult","Microstegium","Microstegium biomass"),each=simtime))
	dfNh$scenario="Invasion"
	
	# combine scenarios
	dfTot=rbind(dfNi,dfNh)
	dfTot$aa=alpha.aa
	dfTot$pp=alpha.pp
	dfTot$ap=alpha.ap
	dfTot$sa=alpha.sa
	
	# return
	return(dfTot)
}

# alpha values
alpha.aa=c(0.1,1) # intraspecific annuals
alpha.pp=c(0.1,1) # intraspecific perennial adults
alpha.ap=c(0.1,1) # p adults on annuals
alpha.sa=c(0.1,1) # annuals on p seedlings
alphas=expand.grid(alpha.aa,alpha.pp,alpha.ap,alpha.sa) #all combinations
colnames(alphas)=c("aa","pp","ap","sa")

# apply function across alpha values
out=list()
for(j in 1:nrow(alphas)){
	out[[j]]=compFun(alphas$aa[j],alphas$sa[j],alphas$ap[j],alphas$pp[j],h.s.comp=F)
}

# dataframe
outdf=data.frame(rbindlist(out))
str(outdf)
outdf$scenario=factor(outdf$scenario,levels=c("Pre-invasion","Invasion"))

# plot
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("PlantDiseaseModel_122018.pdf")
for(k in seq(1,nrow(alphas),by=2)){
	subDat=subset(outdf,pp==alphas$pp[k]&ap==alphas$ap[k]&sa==alphas$sa[k]&species!="Microstegium biomass")
	
	chartTitle=paste("pp:",alphas$pp[k],", ap:",alphas$ap[k],", sa:", alphas$sa[k],", aa:rows",sep="")
	
	print(ggplot(subDat)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_grid(aa~scenario)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")+ggtitle(chartTitle))

}
dev.off()

# repeat simulation with seedling facilitation
out2=list()
for(j in 1:nrow(alphas)){
  out2[[j]]=compFun(alphas$aa[j],alphas$sa[j],alphas$ap[j],alphas$pp[j],h.s.comp=T)
}
# ran with alpha.sp=alpha.pp*10 and seedling population goes to infinity

# dataframe
outdf2=data.frame(rbindlist(out2))
str(outdf2)
outdf2$scenario=factor(outdf2$scenario,levels=c("Pre-invasion","Invasion"))

# plot
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("PlantDiseaseModel_EvSeedlingFac_122018.pdf")
for(k in seq(1,nrow(alphas),by=2)){
  subDat=subset(outdf2,pp==alphas$pp[k]&ap==alphas$ap[k]&sa==alphas$sa[k]&species!="Microstegium biomass")
  
  chartTitle=paste("pp:",alphas$pp[k],", ap:",alphas$ap[k],", sa:", alphas$sa[k],", aa:rows",sep="")
  
  print(ggplot(subDat)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_grid(aa~scenario)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")+ggtitle(chartTitle))
  
}
dev.off()

