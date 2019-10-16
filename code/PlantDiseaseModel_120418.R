# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(data.table)
library(plotly)

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

# seed survival of non-germinated seeds
s.a=1-0.26 #assume 100% viability in first year and a decrease in second year (average value in Huebner 2011, obtained from figure data extraction)

# intraspecific/group competition
alpha.aa=1 # annuals
alpha.ss=0 # perennial seeds on seedling germination/establishment
alpha.pp=1 # perennial adults on seed production

# interspecific/group competition
alpha.as=0 # p seedlings on annuals
alpha.ap=0.1 # p adults on annuals
alpha.sa=0.1 # annuals on p seedling germination/establishment
alpha.sp=alpha.pp*0.01 # p adults on p seedling germination/establishment
alpha.pa=0 # annuals on p adult seed production
alpha.ps=0 # p seedlings on p adult seed production

# adult perennial survival
s.p=1-mean(c(6.2/3,26.3/3,13.2/3))/100 #Malmstrom et al. 2005

# years
simtime=100

# initial numbers
N0.a=1
N0.s=0
N0.p=1


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
	
	Ni.s[t+1]=ifelse(Ni.s[t+1]<1,0,Ni.s[t+1])
	Ni.p[t+1]=ifelse(Ni.p[t+1]<1,0,Ni.p[t+1])
}

# plot
dfNi=data.frame(time=rep(1:simtime,2),N=c(Ni.s,Ni.p),species=rep(c("Elymus seedling","Elymus adult"),each=simtime))
dfNi$scenario="Pre-invasion"
ggplot(dfNi)+geom_line(aes(x=time,y=N,colour=species))+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),legend.position=c(0.8,0.6))+xlab("Years")+ylab("Individuals")


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
	Nh.a[t+1]=s.a*(1-germ.a)*Nh.a[t]+Nh.a[t]*germ.a*lambda.a/(1+alpha.aa*germ.a*Nh.a[t]+alpha.as*germ.s*Nh.s[t]+alpha.ap*Nh.p[t])
	
	Nh.s[t+1]=ifelse(Nh.s[t+1]<1,0,Nh.s[t+1])
	Nh.p[t+1]=ifelse(Nh.p[t+1]<1,0,Nh.p[t+1])
	Nh.a[t+1]=ifelse(Nh.a[t+1]<1,0,Nh.a[t+1])
}

# plot
dfNh=data.frame(time=rep(1:simtime,3),N=c(Nh.s,Nh.p,Nh.a),species=rep(c("Elymus seedling","Elymus adult","Microstegium"),each=simtime))
dfNh$scenario="Invasion"
ggplot(dfNh)+geom_line(aes(x=time,y=N,colour=species))+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),legend.position=c(0.8,0.6))+xlab("Years")+ylab("Individuals")


## Add infection to system

# initialize populations
N.a=rep(NA,simtime)
N.s=rep(NA,simtime)
N.p=rep(NA,simtime)

N.a[1]=N0.a
N.s[1]=Ni.s[simtime]
N.p[1]=Ni.p[simtime]

# infection reduces seed production
lambdaI.a=lambda.a*0.26 # Flory et al. 2011
lambdaI.p=lambda.p*0.9

# assume all individuals get infected

# introduce annual population with infection
for(t in 1:(simtime-1)){
	N.s[t+1]=N.p[t]*lambdaI.p/(1+alpha.pa*germ.a*N.a[t]+alpha.ps*germ.s*N.s[t]+alpha.pp*N.p[t])
	N.p[t+1]=N.s[t]*germ.s/(1+alpha.sa*germ.a*N.a[t]+alpha.ss*germ.s*N.s[t]+alpha.sp*N.p[t])+N.p[t]*s.p
	N.a[t+1]=s.a*(1-germ.a)*N.a[t]+N.a[t]*germ.a*lambdaI.a/(1+alpha.aa*germ.a*N.a[t]+alpha.as*germ.s*N.s[t]+alpha.ap*N.p[t])
	
	N.s[t+1]=ifelse(N.s[t+1]<1,0,N.s[t+1])
	N.p[t+1]=ifelse(N.p[t+1]<1,0,N.p[t+1])
	N.a[t+1]=ifelse(N.a[t+1]<1,0,N.a[t+1])
}

# plot
dfN=data.frame(time=rep(1:simtime,3),N=c(N.s,N.p,N.a),species=rep(c("Elymus seedling","Elymus adult","Microstegium"),each=simtime))
dfN$scenario="Invasion with infection"
ggplot(dfN)+geom_line(aes(x=time,y=N,colour=species))+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),legend.position=c(0.8,0.8))+xlab("Years")+ylab("Individuals")

# plot all three scenarios
dfTot=rbind(dfNi,dfNh,dfN)
dfTot$scenario=factor(dfTot$scenario,levels=c("Pre-invasion","Invasion","Invasion with infection"))
ggplot(dfTot)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_wrap(~scenario)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),legend.position=c(0.55,0.8),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")


## Save N dynamics as a function and run across a range of intraspecific competition coefficients

# simulation function
compFun=function(germ.a,germ.s,lambda.a,lambda.p,s.a,alpha.aa,alpha.ss,alpha.pp,alpha.as,alpha.ap,alpha.sa,alpha.pa,alpha.ps,s.p,simtime,N0.a,N0.s,N0.p,lambdaI.a,lambdaI.p){
	
	# Parameters
	alpha.sp=alpha.pp*0.01
	
	# Perennial starts off by itself

	# initialize populations
	Ni.s=rep(NA,simtime)
	Ni.s[1]=N0.s
	Ni.p=rep(NA,simtime)
	Ni.p[1]=N0.p

	# perennial population dynamics
	for(t in 1:(simtime-1)){
		Ni.s[t+1]=Ni.p[t]*lambda.p/(1+alpha.ps*germ.s*Ni.s[t]+alpha.pp*Ni.p[t])
		Ni.p[t+1]=Ni.s[t]*germ.s/(1+alpha.ss*germ.s*Ni.s[t]+alpha.sp*Ni.p[t])+Ni.p[t]*s.p
		
		Ni.s[t+1]=ifelse(Ni.s[t+1]<1,0,Ni.s[t+1])
		Ni.p[t+1]=ifelse(Ni.p[t+1]<1,0,Ni.p[t+1])
	}

	# save data
	dfNi=data.frame(time=rep(1:simtime,2),N=c(Ni.s,Ni.p),species=rep(c("Elymus seedling","Elymus adult"),each=simtime))
	dfNi$scenario="Pre-invasion"
	
	# Annual invades perennial

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
		Nh.a[t+1]=s.a*(1-germ.a)*Nh.a[t]+Nh.a[t]*germ.a*lambda.a/(1+alpha.aa*germ.a*Nh.a[t]+alpha.as*germ.s*Nh.s[t]+alpha.ap*Nh.p[t])
		
		Nh.s[t+1]=ifelse(Nh.s[t+1]<1,0,Nh.s[t+1])
		Nh.p[t+1]=ifelse(Nh.p[t+1]<1,0,Nh.p[t+1])
		Nh.a[t+1]=ifelse(Nh.a[t+1]<1,0,Nh.a[t+1])
	}

	# save data
	dfNh=data.frame(time=rep(1:simtime,3),N=c(Nh.s,Nh.p,Nh.a),species=rep(c("Elymus seedling","Elymus adult","Microstegium"),each=simtime))
	dfNh$scenario="Invasion"
	
	# Add infection to system

	# initialize populations
	N.a=rep(NA,simtime)
	N.s=rep(NA,simtime)
	N.p=rep(NA,simtime)

	N.a[1]=N0.a
	N.s[1]=Ni.s[simtime]
	N.p[1]=Ni.p[simtime]

	# introduce annual population with infection
	for(t in 1:(simtime-1)){
		N.s[t+1]=N.p[t]*lambdaI.p/(1+alpha.pa*germ.a*N.a[t]+alpha.ps*germ.s*N.s[t]+alpha.pp*N.p[t])
		N.p[t+1]=N.s[t]*germ.s/(1+alpha.sa*germ.a*N.a[t]+alpha.ss*germ.s*N.s[t]+alpha.sp*N.p[t])+N.p[t]*s.p
		N.a[t+1]=s.a*(1-germ.a)*N.a[t]+N.a[t]*germ.a*lambdaI.a/(1+alpha.aa*germ.a*N.a[t]+alpha.as*germ.s*N.s[t]+alpha.ap*N.p[t])
		
		N.s[t+1]=ifelse(N.s[t+1]<1,0,N.s[t+1])
		N.p[t+1]=ifelse(N.p[t+1]<1,0,N.p[t+1])
		N.a[t+1]=ifelse(N.a[t+1]<1,0,N.a[t+1])
	}

	# save data
	dfN=data.frame(time=rep(1:simtime,3),N=c(N.s,N.p,N.a),species=rep(c("Elymus seedling","Elymus adult","Microstegium"),each=simtime))
	dfN$scenario="Invasion with infection"
	
	# combine all three scenarios
	dfTot=rbind(dfNi,dfNh,dfN)
	dfTot$aa=alpha.aa
	dfTot$pp=alpha.pp
	dfTot$ap=alpha.ap
	dfTot$sa=alpha.sa
	
	# return
	return(dfTot)
}

# alpha values
alpha.aa=c(0.1,1,2) # intraspecific annuals
alpha.pp=c(0.1,1,2) # intraspecific perennial adults
alpha.ap=c(0.1,1,2) # p adults on annuals
alpha.sa=c(0.1,1,2) # annuals on p seedlings
alphas=expand.grid(alpha.aa,alpha.pp,alpha.ap,alpha.sa) #all combinations
colnames(alphas)=c("aa","pp","ap","sa")

# apply function across alpha values
out=list()
for(j in 1:nrow(alphas)){
	out[[j]]=compFun(germ.a,germ.s,lambda.a,lambda.p,s.a,alphas$aa[j],alpha.ss,alphas$pp[j],alpha.as,alphas$ap[j],alphas$sa[j],alpha.pa,alpha.ps,s.p,simtime,N0.a,N0.s,N0.p,lambdaI.a,lambdaI.p)
	
}

# dataframe
outdf=data.frame(rbindlist(out))
str(outdf)
outdf$scenario=factor(outdf$scenario,levels=c("Pre-invasion","Invasion","Invasion with infection"))

# plot
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("PlantDiseaseModel_120418.pdf")
for(k in seq(1,nrow(alphas),by=3)){
	subDat=subset(outdf,pp==alphas$pp[k]&ap==alphas$ap[k]&sa==alphas$sa[k])
	
	chartTitle=paste("pp:",alphas$pp[k],", ap:",alphas$ap[k],", sa:", alphas$sa[k],", aa:rows",sep="")
	
	print(ggplot(subDat)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_grid(aa~scenario)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")+ggtitle(chartTitle))

}
dev.off()

# final time point
findf=subset(outdf,time==simtime)

# when does invasion occur in the absence of infection?
invdf=subset(findf,species=="Microstegium"&N>0)
dim(invdf)
invdf #ap=0.1 and pp>0.1
#none of these have infection, so infection reduces invasion every time

# plot
plot_ly(invdf,x=~aa,y=~pp,z=~N,color=~sa)%>%
add_markers()%>%
layout(scene=list(xaxis=list(title="annual intraspecific"),yaxis=list(title="perennial intraspecific"),zaxis=list(title="Mv density")))
#annual intraspecific is dribing most of the variation

# plot
ggplot(invdf,aes(x=aa,y=N,colour=as.factor(pp)))+geom_point()+geom_line()+facet_wrap(~sa)
