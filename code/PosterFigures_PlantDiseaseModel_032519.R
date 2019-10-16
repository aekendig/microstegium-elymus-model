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
axisText=20
axisTitle=23
legendText=20
legendTitle=23

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
lambda.a=512
lambda.p=72
lambda.s=9

# seed production competition parameters (greenhouse/lit)
alpha.aa=mean(c(0.015,0.001))
alpha.as=0.0037
alpha.sa=0.091
alpha.ps=0
alpha.pp=0.06

alpha.ap=alpha.as*10
alpha.ss=alpha.as
alpha.sp=alpha.ap
alpha.pa=alpha.sa/10

# convert seed production to biomass
c.a=30/6500

# decrease seed production due to infection
tolMin=0.31
tol.a=0.6
tol.p=0.95

# years
simtime.p=100
simtime.i=300

# initial numbers
N0.a=100
N0.s=0
N0.p=1

## Save N dynamics as a function and run across a range of intraspecific competition coefficients

# simulation function
compFun=function(lambda.a, m.p, tol.a, tol.p, h.s.fac){
	
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
	  h.s=ifelse(h.s.fac==T,h.s_fun(Ni.s[t]),h.s.mean)
	  h.s=ifelse(is.na(h.s)&h.s.fac==T,1,h.s)
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
	N.a=rep(NA,simtime.i)
	N.s=rep(NA,simtime.i)
	N.p=rep(NA,simtime.i)
	L=rep(NA,simtime.i)

	N.a[1]=N0.a
	N.s[1]=Ni.s[simtime.p]
	N.p[1]=Ni.p[simtime.p]
	L[1]=0

	# introduce annual population
	for(t in 1:(simtime.i-1)){	
		
		# infection-dependent seed production
		lam.a=ifelse(t>100,tol.a*lambda.a,lambda.a)
		lam.p=ifelse(t>100,tol.p*lambda.p,lambda.p)
		lam.s=ifelse(t>100,tol.p*lambda.s,lambda.s)	
		
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

	# save data
	dfN=data.frame(time=rep(1:simtime.i,4),N=c(N.s,N.p,N.a,L),species=rep(c("Elymus seedling","Elymus adult","Microstegium","Microstegium litter"),each=simtime.i))
	dfN$scenario="Invasion"

	
	# combine scenarios
	dfTot=rbind(dfNi,dfN)
	dfTot$lambda.a=lambda.a
	dfTot$m.p=m.p
	dfTot$tol.a=tol.a
	dfTot$tol.p=tol.p
	
	# return
	return(dfTot)
}


## Examine model outcome with default values

outDef = compFun(lambda.a = lambda.a, m.p = m.p, tol.a = tol.a, tol.p = tol.p, h.s.fac = F)
outDef$scenario=factor(outDef$scenario,levels=c("Pre-invasion","Invasion"))

# plot
ggplot(outDef)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_wrap(~scenario)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")
# invasion quickly suppressed by disease

# increase propagule pressure

outHi = compFun(lambda.a = 10000, m.p = m.p, tol.a = tol.a, tol.p = tol.p, h.s.fac = F)
outHi$scenario=factor(outHi$scenario,levels=c("Pre-invasion","Invasion"))

# plot
ggplot(outHi)+geom_line(aes(x=time,y=N,colour=species),size=1)+facet_wrap(~scenario)+theme(axis.text=element_text(size=12,colour="black"),axis.title=element_text(size=14,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=12),legend.title=element_text(size=12),legend.key=element_blank(),strip.background=element_blank(),strip.text=element_text(size=12))+xlab("Years")+ylab("Individuals")
# invasion worked


## Examine model outcomes with a range of community composition values

# perennial parameters
annFec=c(1,10,100,500,1000,5000,10000,50000,1e5,5e5,1e6) # initial annual densities
perSur=seq(0.7,1,length.out=10) # perennial survival
params_LaMp=expand.grid(annFec, perSur) #all combinations
colnames(params_LaMp)=c("lambda.a","m.p")

# repeat simulation over parameter values
outSim_LaMp=list()
for(j in 1:nrow(params_LaMp)){
  outSim_LaMp[[j]]=compFun(lambda.a=params_LaMp$lambda.a[j],m.p=params_LaMp$m.p[j],tol.a=tol.a,tol.p=tol.p,h.s.fac=F)
}

# dataframe
outDF_LaMp=data.frame(rbindlist(outSim_LaMp))
str(outDF_LaMp)
outDF_LaMp$scenario=factor(outDF_LaMp$scenario,levels=c("Pre-invasion","Invasion"))

# final day dataframe
finDF_LaMp=subset(outDF_LaMp,time==simtime.i)%>%spread(.,species,N)%>%rename(.,EvA="Elymus adult",EvS="Elymus seedling",Mv="Microstegium",MvL="Microstegium litter")
finDF_LaMp$propEv=with(finDF_LaMp,(EvA+EvS)/(EvA+EvS+Mv))
finDF_LaMp$propMv=with(finDF_LaMp,Mv/(EvA+EvS+Mv))

# plot
ggplot(finDF_LaMp,aes(x=log10(lambda.a),y=m.p,colour=propEv))+geom_point(size=5,shape=15)+scale_colour_gradientn(colors=c("#92D050","white","#002060"))

## Examine model outcomes with a range of tolerance values

# perennial parameters
annTol=seq(0.01,1,length.out=10)
perTol=seq(0.01,1,length.out=10)
params_tol=expand.grid(annTol, perTol) #all combinations
colnames(params_tol)=c("tol.a","tol.p")

# repeat simulation over parameter values
outSim_tol=list()
for(j in 1:nrow(params_tol)){
  outSim_tol[[j]]=compFun(lambda.a=lambda.a,m.p=m.p,tol.a=params_tol$tol.a[j],tol.p=params_tol$tol.p[j],h.s.fac=F)
}

# dataframe
outDF_tol=data.frame(rbindlist(outSim_tol))
str(outDF_tol)
outDF_tol$scenario=factor(outDF_tol$scenario,levels=c("Pre-invasion","Invasion"))

# final day dataframe
finDF_tol=subset(outDF_tol,time==simtime.i)%>%spread(.,species,N)%>%rename(.,EvA="Elymus adult",EvS="Elymus seedling",Mv="Microstegium",MvL="Microstegium litter")
finDF_tol$propEv=with(finDF_tol,(EvA+EvS)/(EvA+EvS+Mv))
finDF_tol$propMv=with(finDF_tol,Mv/(EvA+EvS+Mv))

# plot
ggplot(finDF_tol,aes(x=tol.a,y=tol.p,colour=propEv))+geom_point(size=5,shape=15)+scale_colour_gradientn(colors=c("#92D050","white","#002060"))


## Final figure

# scale by max N
maxN=max(outDF_tol$N)

# choose specific simulations
aexcl=subset(outDF_tol,tol.p==0.01&tol.a==1&scenario=="Invasion"&species!="Microstegium litter")%>%spread(.,species,N)%>%rename(.,EvA="Elymus adult",EvS="Elymus seedling",Mv="Microstegium")%>%mutate(Ev=EvA+EvS)%>%gather(key=species,value=N,-c(time:tol.p))%>%filter(species%in%c("Ev","Mv"))
aexcl$Nscal=aexcl$N/maxN
max(aexcl$Nscal) #1

pexcl=subset(outDF_tol,tol.p==1&tol.a==0.01&scenario=="Invasion"&species!="Microstegium litter")%>%spread(.,species,N)%>%rename(.,EvA="Elymus adult",EvS="Elymus seedling",Mv="Microstegium")%>%mutate(Ev=EvA+EvS)%>%gather(key=species,value=N,-c(time:tol.p))%>%filter(species%in%c("Ev","Mv"))
pexcl$Nscal=pexcl$N/maxN
max(pexcl$Nscal) #0.06

coex=subset(outDF_tol,tol.p==0.45&tol.a==0.67&scenario=="Invasion"&species!="Microstegium litter")%>%spread(.,species,N)%>%rename(.,EvA="Elymus adult",EvS="Elymus seedling",Mv="Microstegium")%>%mutate(Ev=EvA+EvS)%>%gather(key=species,value=N,-c(time:tol.p))%>%filter(species%in%c("Ev","Mv"))
coex$Nscal=coex$N/maxN
max(coex$Nscal) #0.06
coex$species=revalue(as.factor(coex$species),c("Ev"="native","Mv"="invader"))

plot1=ggplot(aexcl,aes(x=time,y=Nscal,colour=species))+geom_line(size=2)+ylim(0,1)+scale_colour_manual(values=c("#407879","black"),name="Population")+theme(axis.text=element_text(size=axisText,colour="black"),axis.title=element_text(size=axisTitle,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=legendText),legend.title=element_text(size=legendTitle),legend.key=element_blank(),legend.position="none",plot.title=element_text(size=axisTitle))+xlab("Years")+ylab("Scaled population size")+ggtitle("Invader wins\n(higher disease tolerance)")

plot2=ggplot(coex,aes(x=time,y=Nscal,colour=species))+geom_line(size=2)+ylim(0,0.06)+scale_colour_manual(values=c("#407879","black"),name="Species")+theme(axis.text=element_text(size=axisText,colour="black"),axis.title=element_text(size=axisTitle,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=legendText),legend.title=element_text(size=legendTitle),legend.key=element_blank(),legend.position=c(0.4,0.85),plot.title=element_text(size=axisTitle))+xlab("Years")+ylab("Scaled population size")+ggtitle("Coexistence\n(similar disease tolerance)")

plot3=ggplot(pexcl,aes(x=time,y=Nscal,colour=species))+geom_line(size=2)+ylim(0,0.06)+scale_colour_manual(values=c("#407879","black"),name="Population")+theme(axis.text=element_text(size=axisText,colour="black"),axis.title=element_text(size=axisTitle,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=legendText),legend.title=element_text(size=legendTitle),legend.key=element_blank(),legend.position="none",plot.title=element_text(size=axisTitle))+xlab("Years")+ylab("Scaled population size")+ggtitle("Native wins\n(higher disease tolerance)")

plot4=ggplot(finDF_tol,aes(x=tol.a,y=tol.p,colour=propEv))+geom_point(size=30,shape=15)+scale_colour_gradientn(colors=c("black","white","#407879"),name="Native\nrelative\nabundance")+theme(axis.text=element_text(size=axisText,colour="black"),axis.title=element_text(size=axisTitle,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=legendText),legend.title=element_text(size=legendTitle),legend.key=element_blank(),legend.position="right")+xlab("Invader disease tolerance")+ylab("Native disease tolerance")

plot5=ggplot(finDF_LaMp,aes(x=log10(lambda.a),y=m.p,colour=propEv))+geom_point(size=45,shape=15)+scale_colour_gradientn(colors=c("black","white","#407879"))+theme(axis.text=element_text(size=axisText,colour="black"),axis.title=element_text(size=axisTitle,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA),legend.text=element_text(size=legendText),legend.title=element_text(size=legendTitle),legend.key=element_blank(),legend.position="none")+xlab(expression(paste(log[10],"(Invader seed prod.)",sep="")))+ylab("Native adult survival")


setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("CompetitionOutcomes_PlantDiseaseModel_032619.pdf",width=25,height=5)
plot_grid(plot1,plot2,plot3,plot4,plot5,nrow=1,rel_widths=c(1,1,1,1.45,1))
dev.off()

