## need to update model to match model building powerpoint

# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(data.table)
library(plotly)
library(tidyr)
library(viridis)

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
alpha.aa1=0.05
alpha.aa2=0.15
alpha.as=0.0037
alpha.ap=0.037
alpha.sa=0.091
alpha.ss=0.0037
alpha.sp=0.037
alpha.pa=0.0091
alpha.ps=0
alpha.pp1=0.02
alpha.pp2=0.06

# convert seed production to biomass
c.a=30/6500

# decrease seed production due to infection
virMax=0.31

# years
simtime.p=100
simtime.i=300

# initial numbers
N0.a=1
N0.s=0
N0.p=1

## Save N dynamics as a function and run across a range of intraspecific competition coefficients

# simulation function
compFun=function(alpha.aa,alpha.pp,h.s.comp,vir.a,vir.p){
	
	# functions
  h.s_fun=function(x){
    h.out=exp(h.s.int+h.s.slope*x)/(1+exp(h.s.int+h.s.slope*x))
    return(h.out)
  }
  
	# perennial starts off by itself

	# initialize populations
	Ni.s=rep(NA,simtime.p)
	Ni.s[1]=N0.s
	Ni.p=rep(NA,simtime.p)
	Ni.p[1]=N0.p

	# perennial population dynamics
	for(t in 1:(simtime.p-1)){
	  
	  # composite parameters
	  h.s=ifelse(h.s.comp==T,h.s_fun(Ni.s[t]),h.s.mean)
	  h.s=ifelse(is.na(h.s)&h.s.comp==T,1,h.s)
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
	dfNi=data.frame(time=rep(1:simtime.p,2),N=c(Ni.s,Ni.p),species=rep(c("Elymus seedling","Elymus adult"),each=simtime.p))
	dfNi$scenario="Pre-invasion"

	# annual invades perennial

	# initialize populations
	Nh.a=rep(NA,simtime.i)
	Nh.s=rep(NA,simtime.i)
	Nh.p=rep(NA,simtime.i)
	Bh.a=rep(NA,simtime.i)

	Nh.a[1]=N0.a
	Nh.s[1]=Ni.s[simtime.p]
	Nh.p[1]=Ni.p[simtime.p]
	Bh.a[1]=0

	# introduce annual population
	for(t in 1:(simtime.i-1)){	
		
		# infection-dependent seed production
		lam.a=ifelse(t>100,vir.a*lambda.a,lambda.a)
		lam.p=ifelse(t>100,vir.p*lambda.p,lambda.p)
		lam.s=ifelse(t>100,vir.p*lambda.s,lambda.s)	
		
		# composite paramters
	  	h.s=ifelse(h.s.comp==T,h.s_fun(Ni.s[t]),h.s.mean)
	  	h.s=ifelse(is.na(h.s)&h.s.comp==T,1,h.s)
	  	Lh=b*Bh.a[t]
	  	g.s=gamma.s/(1+alpha.sL*Lh)
	  	g.a=gamma.a/(1+alpha.aL*Lh)
	  	f.s=lam.s/(1+alpha.ss*g.s*h.s*Nh.s[t]+alpha.sp*m.p*Nh.p[t]+alpha.sa*g.a*h.a*Nh.a[t])
	  	f.p=lam.p/(1+alpha.ps*g.s*h.s*Nh.s[t]+alpha.pp*m.p*Nh.p[t]+alpha.pa*g.a*h.a*Nh.a[t])
	  	f.a=lam.a/(1+alpha.as*g.s*h.s*Nh.s[t]+alpha.ap*m.p*Nh.p[t]+alpha.aa*g.a*h.a*Nh.a[t])

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


## Examine model outcomes with a few alpha values

# alpha values
alpha.aa=c(alpha.aa1,alpha.aa2) # intraspecific annuals
alpha.pp=c(alpha.pp1,alpha.pp2) # intraspecific perennial adults
alphas=expand.grid(alpha.aa,alpha.pp) #all combinations
colnames(alphas)=c("aa","pp")

# apply function across alpha values
out=list()
for(j in 1:nrow(alphas)){
	out[[j]]=compFun(alphas$aa[j],alphas$pp[j],h.s.comp=F,vir.a=1,vir.p=1)
}

# dataframe
outdf=data.frame(rbindlist(out))
str(outdf)
outdf$scenario=factor(outdf$scenario,levels=c("Pre-invasion","Invasion"))

# plot
ggplot(subset(outdf,scenario=="Invasion"))+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_grid(aa~pp)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")

# repeat simulation with seedling facilitation
out2=list()
for(j in 1:nrow(alphas)){
  out2[[j]]=compFun(alphas$aa[j],alphas$pp[j],h.s.comp=T,vir.a=1,vir.p=1)
}

# dataframe
outdf2=data.frame(rbindlist(out2))
str(outdf2)
outdf2$scenario=factor(outdf2$scenario,levels=c("Pre-invasion","Invasion"))

# plot
ggplot(subset(outdf2,scenario=="Invasion"))+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_grid(aa~pp)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")
# Ev facilitation can prevent invasion

# repeat simulation with Mv infection
out3=list()
for(j in 1:nrow(alphas)){
  out3[[j]]=compFun(alphas$aa[j],alphas$pp[j],h.s.comp=F,vir.a=virMax,vir.p=1)
}

# dataframe
outdf3=data.frame(rbindlist(out3))
str(outdf3)
outdf3$scenario=factor(outdf3$scenario,levels=c("Pre-invasion","Invasion"))

# plot
ggplot(subset(outdf3,scenario=="Invasion"))+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_grid(aa~pp)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")
# Mv population crashes

# repeat simulation with Mv and Ev infection
out4=list()
for(j in 1:nrow(alphas)){
  out4[[j]]=compFun(alphas$aa[j],alphas$pp[j],h.s.comp=F,vir.a=virMax,vir.p=virMax)
}

# dataframe
outdf4=data.frame(rbindlist(out4))
str(outdf4)
outdf4$scenario=factor(outdf4$scenario,levels=c("Pre-invasion","Invasion"))

# plot
ggplot(subset(outdf4,scenario=="Invasion"))+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_grid(aa~pp)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")
# Mv and Ev coexist


## Examine model outcomes with a range of intraspecific competition and virulence for Ev

# perennial parameters
alpha.pp2=seq(0.001,1,length.out=10) # intraspecific perennial adults
virulence.p=seq(0.1,1,length.out=10) # reduction in seed production due to infection for perennials
params.p=expand.grid(alpha.pp2,virulence.p) #all combinations
colnames(params.p)=c("pp","vir")

# constant annual parameters
alpha.aa=mean(c(0.015,0.001))
vir.a=0.6

# repeat simulation over perennial parameters
out5=list()
for(j in 1:nrow(params.p)){
  out5[[j]]=compFun(alpha.aa=alpha.aa,alpha.pp=params.p$pp[j],h.s.comp=F,vir.a=vir.a,vir.p=params.p$vir[j])
}

# dataframe
outdf5=data.frame(rbindlist(out5))
str(outdf5)
outdf5$scenario=factor(outdf5$scenario,levels=c("Pre-invasion","Invasion"))

# print plots to check 
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("PlantDiseaseModel_PerennialParamSims_020519.pdf")
for(k in 1:nrow(params.p)){
  subDat=subset(outdf5,pp==params.p$pp[k]&vir.p==params.p$vir[k])
  
  chartTitle=paste("pp:",params.p$pp[k],", vir.p:",params.p$vir[k],sep="")
  
  print(ggplot(subDat)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_wrap(~scenario,nrow=1,scales="free_x")+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")+ggtitle(chartTitle)
  )
  
}
dev.off()


## Repeat above with finer scale of parameters

# perennial parameters
alpha.pp3=seq(0.001,0.3,length.out=40) # intraspecific perennial adults
virulence.p2=seq(0.1,1,length.out=40) # reduction in seed production due to infection for perennials
params.p2=expand.grid(alpha.pp3,virulence.p2) #all combinations
colnames(params.p2)=c("pp","vir")

# repeat simulation over perennial parameters
out6=list()
for(j in 1:nrow(params.p2)){
  out6[[j]]=compFun(alpha.aa=alpha.aa,alpha.pp=params.p2$pp[j],h.s.comp=F,vir.a=vir.a,vir.p=params.p2$vir[j])
}

# dataframe
outdf6=data.frame(rbindlist(out6))
str(outdf6)
outdf6$scenario=factor(outdf6$scenario,levels=c("Pre-invasion","Invasion"))

# final day dataframe
findf6=subset(outdf6,time==simtime.i)%>%spread(.,species,N)%>%rename(.,EvA="Elymus adult",EvS="Elymus seedling",Mv="Microstegium",MvB="Microstegium biomass")
findf6$propEv=with(findf6,(EvA+EvS)/(EvA+EvS+Mv))
findf6$propMv=with(findf6,Mv/(EvA+EvS+Mv))

# plot
ggplot(findf6,aes(x=pp,y=vir.p,colour=propEv))+geom_point(size=5,shape=15)+scale_colour_gradientn(colors=c("#92D050","white","#002060"))


## Repeat again with higher virulence for Mv

# repeat simulation over perennial parameters
out7=list()
for(j in 1:nrow(params.p2)){
  out7[[j]]=compFun(alpha.aa=alpha.aa,alpha.pp=params.p2$pp[j],h.s.comp=F,vir.a=virMax,vir.p=params.p2$vir[j])
}

# dataframe
outdf7=data.frame(rbindlist(out7))
str(outdf7)
outdf7$scenario=factor(outdf7$scenario,levels=c("Pre-invasion","Invasion"))

# final day dataframe
findf7=subset(outdf7,time==simtime.i)%>%spread(.,species,N)%>%rename(.,EvA="Elymus adult",EvS="Elymus seedling",Mv="Microstegium",MvB="Microstegium biomass")
findf7$propEv=with(findf7,(EvA+EvS)/(EvA+EvS+Mv))
findf7$propMv=with(findf7,Mv/(EvA+EvS+Mv))

# plot
ggplot(findf7,aes(x=pp,y=vir.p,colour=propEv))+geom_point(size=5,shape=15)+scale_colour_gradientn(colors=c("#92D050","white","#002060"))


## Final figure

# scale by max N
maxN=max(outdf6$N)

aexcl=subset(outdf6,round(vir.p,2)==0.61&round(pp,4)==0.0087&scenario=="Invasion"&species!="Microstegium biomass")
aexcl$Nscal=aexcl$N/maxN
max(aexcl$Nscal) #0.12

pexcl=subset(outdf6,round(vir.p,2)==0.61&round(pp,4)==0.3000&scenario=="Invasion"&species!="Microstegium biomass")
pexcl$Nscal=pexcl$N/maxN
max(pexcl$Nscal) #0.11

coex=subset(outdf6,round(vir.p,2)==0.61&round(pp,4)==0.0777&scenario=="Invasion"&species!="Microstegium biomass")
coex$Nscal=coex$N/maxN
max(coex$Nscal) #0.016

ggplot(aexcl,aes(x=time,y=Nscal,colour=species))+geom_line(size=1)+ylim(0,0.12)

ggplot(pexcl,aes(x=time,y=Nscal,colour=species))+geom_line(size=1)+ylim(0,0.12)

ggplot(coex,aes(x=time,y=Nscal,colour=species))+geom_line(size=1)+ylim(0,0.12)