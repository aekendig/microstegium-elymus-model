# clear all existing data
rm(list=ls())

# open libraries
library(data.table)
library(plotly)
library(cowplot)
library(popbio)
library(tidyverse)

# plotting parameters
axisText=12
axisTitle=14
legendText=12
legendTitle=0

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
s.a=0.74
s.s=0.76

# germination rate
gamma.a=0.7
gamma.s=0.66

# germination reduction with litter (per gram)
alpha.aL=alpha.sL=-0.35/374 #-0.0009

# proportion of Mv biomass that becomes litter the next year
b=0.56

# survival to reproduction
h.a=0.95
h.s=0.40
h.p=0.83

# seed production rates
lambda.a=6500 
lambda.p=435

# seed production competition parameters
#alpha.aa=mean(c(0.015,0.001)) #0.008
alpha.aa=0.015
#alpha.sa=mean(c(0.724,9.574)) #5.149 - Mv can't establish invasion
alpha.sa=9.574
alpha.pa=alpha.sa/10 #0.9574

#alpha.as=mean(c(0.054,0.357)) #0.2055 - Mv can't establish invasion
alpha.as=0
alpha.ss=mean(c(0.002,0.006)) #0.004
alpha.ps=0

alpha.ap=mean(c(0.054,0.357)) #0.2055
alpha.sp=alpha.ap #0.2055
alpha.pp=alpha.sp*10 #2.055
	
	
# convert seed production to biomass
c.a=30/6500

# decrease seed production due to infection
tolMax=0.31
tolMin=0.6

# parameter list
params=c(m.p,s.a,s.s,gamma.a,gamma.s,alpha.aL,alpha.sL,b,h.a,h.s,h.p,lambda.a,lambda.p,alpha.aa,alpha.as,alpha.ap,alpha.sa,alpha.ss,alpha.sp,alpha.pa,alpha.ps,alpha.pp,c.a)

# parameter names
paramNames=c("m.p","s.a","s.s","gamma.a","gamma.s","alpha.aL","alpha.sL","b","h.a","h.s","h.p","lambda.a","lambda.p","alpha.aa","alpha.as","alpha.ap","alpha.sa","alpha.ss","alpha.sp","alpha.pa","alpha.ps","alpha.pp","c.a","tol.a","tol.p")


## Population Model

simFun=function(params,invader,N0.a,N0.s,N0.p,L0,simtime){
	
	# define parameters
	m.p=params[1]
	s.a=params[2]
	s.s=params[3]
	gamma.a=params[4]
	gamma.s=params[5]
	alpha.aL=params[6]
	alpha.sL=params[7]
	b=params[8]	
	h.a=params[9]
	h.s=params[10]
	h.p=params[11]
	lambda.a=params[12]
	lambda.p=params[13]
	alpha.aa=params[14]
	alpha.as=params[15]
	alpha.ap=params[16]
	alpha.sa=params[17]
	alpha.ss=params[18]
	alpha.sp=params[19]
	alpha.pa=params[20]
	alpha.ps=params[21]
	alpha.pp=params[22]
	c.a=params[23]
	tol.a=params[24]
	tol.p=params[25]
	
	# calculate parameters
	lambda.s=lambda.p/10
	lam.a=tol.a*lambda.a
	lam.p=tol.p*lambda.p
	lam.s=tol.p*lambda.s	

	# initialize populations
	N.a=rep(NA,simtime)
	N.s=rep(NA,simtime)
	N.p=rep(NA,simtime)
	L=rep(NA,simtime)

	N.a[1]=N0.a
	N.s[1]=N0.s
	N.p[1]=N0.p
	L[1]=L0
	
	# initialize grwr
	grwr=rep(NA,simtime)

	# simulate population dynamics
	for(t in 1:(simtime-1)){	
		
		# composite paramters
	  	g.s=gamma.s+(alpha.sL*L[t])
	  	g.s=ifelse(g.s<0, 0, g.s)
	  	g.a=gamma.a+(alpha.aL*L[t])
	  	g.a=ifelse(g.a<0, 0, g.a)
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
		
		# perennial transition matrix
		pmat=matrix(c(s.s*(1-g.s)+g.s*h.s*f.s,m.p*f.p,g.s*h.s,m.p),nrow=2,byrow=T)
		grwr.p=eigen.analysis(pmat)$lambda1
		
		# annual growth rate when rare
		grwr.a=s.a*(1-g.a)+g.a*h.a*f.a
		
		# save grwr
		grwr[t]=ifelse(invader=="P",grwr.p,grwr.a)
	}
	
	# final grwr
	
	# perennial transition matrix
	t = simtime
	pmat=matrix(c(s.s*(1-g.s)+g.s*h.s*f.s,m.p*f.p,g.s*h.s,m.p),nrow=2,byrow=T)
	grwr.p=eigen.analysis(pmat)$lambda1
		
	# annual growth rate when rare
	grwr.a=s.a*(1-g.a)+g.a*h.a*f.a
		
	# save grwr
	grwr[t]=ifelse(invader=="P",grwr.p,grwr.a)
	
	# save data
	dfN=data.frame(time=rep(1:simtime,4),N=c(N.s,N.p,N.a,L),species=rep(c("Elymus seedling","Elymus adult","Microstegium","Microstegium litter"),each=simtime))
	dfN$grwr=grwr
	
	# return
	return(dfN)
}


## Simulations

# simulation time
years = 500

# Annual alone, no infection
a_alone_none=simFun(params=c(params, 1, 1),invader="A",N0.a=1,N0.s=0,N0.p=0,L0=0,simtime=years)

ggplot(a_alone_none, aes(x = time, y = N, color = species)) +
geom_line()

filter(a_alone_none, time == years)

# Annual alone, high tolerance
a_alone_high=simFun(params=c(params, tolMax, tolMax),invader="A",N0.a=1,N0.s=0,N0.p=0,L0=0,simtime=years)

ggplot(a_alone_high, aes(x = time, y = N, color = species)) +
geom_line()

filter(a_alone_high, time == years)

# Perennial alone, no infection
p_alone_none=simFun(params=c(params, 1, 1),invader="P",N0.a=0,N0.s=1,N0.p=0,L0=0,simtime=years)

p_alone_none_2 <- p_alone_none %>%
	spread(.,species,N) %>% 
	rename(.,EvA="Elymus adult",EvS="Elymus seedling")%>% 
	mutate(Ev=EvA+EvS)%>% 
	gather(key=species,value=N,-c(time,grwr)) %>%
	mutate(species = recode(species, "Ev" = "Native"))

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("NativeAlone_PlantDiseaseModel_080519.pdf",width=3,height=3)
p_alone_none_2 %>%
	filter(species == "Native" & time <= 200) %>%
	ggplot(aes(x = time, y = N, color = species)) +
	geom_line(size = 2)+
	scale_colour_manual(values=c("#407879"), name = "") +
	theme_bw() +
	theme(axis.text=element_text(size=axisText,colour="black"),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position=c(0.7, 0.4),
	plot.title=element_text(size=axisTitle)) +
	xlab("Years") +
	ylab("Population size")
dev.off()

filter(p_alone_none, time == years)
filter(p_alone_none, time == 200)

# Annual invades perennial, no infection
N0.s.pa = filter(p_alone_none, time == years & species == "Elymus seedling")$N
N0.p.pa = filter(p_alone_none, time == years & species == "Elymus adult")$N

a_inv_none=simFun(params=c(params, 1, 1),invader="A",N0.a=1,N0.s=N0.s.pa,N0.p=N0.p.pa,L0=0,simtime=years)

a_inv_none_2 <- a_inv_none %>%
	spread(.,species,N) %>% 
	rename(.,EvA="Elymus adult",EvS="Elymus seedling")%>% 
	mutate(Ev=EvA+EvS)%>% 
	gather(key=species,value=N,-c(time,grwr)) %>%
	mutate(species = recode(species, "Ev" = "Native", "Microstegium" = "Invader", "Microstegium litter" = "Invader litter"),
	time = time + 200)

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("AnnualInvades_PlantDiseaseModel_080519.pdf",width=3,height=3)
a_inv_none_2 %>%
	filter(!(species %in% c("EvS", "EvA")) & time <= 400) %>%
	ggplot(aes(x = time, y = N, color = species)) +
	geom_line(size = 1)+
	scale_colour_manual(values=c("black", "gray", "#407879"), name = "") +
	theme_bw() +
	theme(axis.text=element_text(size=axisText,colour="black"),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position=c(0.6, 0.3),
	plot.title=element_text(size=axisTitle)) +
	xlab("Years") +
	ylab("Population size")
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("AnnualInvades_NoAnn_PlantDiseaseModel_080519.pdf",width=3,height=3)
a_inv_none_2 %>%
	filter(!(species %in% c("EvS", "EvA", "Invader")) & time <= 400) %>%
	ggplot(aes(x = time, y = N, color = species)) +
	geom_line(size = 1)+
	scale_colour_manual(values=c("gray", "#407879"), name = "") +
	theme_bw() +
	theme(axis.text=element_text(size=axisText,colour="black"),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position="none",
	plot.title=element_text(size=axisTitle)) +
	xlab("Years") +
	ylab("Population size")
dev.off()

filter(a_inv_none, time == years)

# Infection spreads, annual has low tolerance
N0.a.ai = mean(filter(a_inv_none, species == "Microstegium" & time > 400)$N)
L0.ai = mean(filter(a_inv_none, species == "Microstegium litter" & time > 400)$N)
N0.s.ai = mean(filter(a_inv_none, species == "Elymus seedling" & time > 400)$N)
N0.p.ai = mean(filter(a_inv_none, species == "Elymus adult" & time > 400)$N)

inf_low_tol=simFun(params=c(params, tolMin, 1),invader="A",N0.a=N0.a.ai,N0.s=N0.s.ai,N0.p=N0.p.ai,L0=L0.ai,simtime=years)

inf_low_tol_2 <- inf_low_tol %>%
	spread(.,species,N) %>% 
	rename(.,EvA="Elymus adult",EvS="Elymus seedling")%>% 
	mutate(Ev=EvA+EvS)%>% 
	gather(key=species,value=N,-c(time,grwr)) %>%
	mutate(species = recode(species, "Ev" = "Native", "Microstegium" = "Invader", "Microstegium litter" = "Invader litter"),
	time = time + 400)

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("AnnualInfection_PlantDiseaseModel_080519.pdf",width=3,height=3)
inf_low_tol_2 %>%
	filter(!(species %in% c("EvS", "EvA")) & time <= 600) %>%
	ggplot(aes(x = time, y = N, color = species)) +
	geom_line(size = 0.5)+
	scale_colour_manual(values=c("black", "gray", "#407879"), name = "") +
	theme_bw() +
	theme(axis.text=element_text(size=axisText,colour="black"),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position="none",
	plot.title=element_text(size=axisTitle)) +
	xlab("Years") +
	ylab("Population size")
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("AnnualInfection_NoAnn_PlantDiseaseModel_080519.pdf",width=3,height=3)
inf_low_tol_2 %>%
	filter(!(species %in% c("EvS", "EvA", "Invader")) & time <= 600) %>%
	ggplot(aes(x = time, y = N, color = species)) +
	geom_line(size = 1)+
	scale_colour_manual(values=c("gray", "#407879"), name = "") +
	theme_bw() +
	theme(axis.text=element_text(size=axisText,colour="black"),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position="none",
	plot.title=element_text(size=axisTitle)) +
	xlab("Years") +
	ylab("Population size")
dev.off()

filter(inf_low_tol, time == years)

# Infection spreads, both have low tolerance
inf_low_tol_both=simFun(params=c(params, tolMin, tolMin),invader="A",N0.a=N0.a.ai,N0.s=N0.s.ai,N0.p=N0.p.ai,L0=L0.ai,simtime=years)

inf_low_tol_both_2 <- inf_low_tol_both %>%
	spread(.,species,N) %>% 
	rename(.,EvA="Elymus adult",EvS="Elymus seedling")%>% 
	mutate(Ev=EvA+EvS)%>% 
	gather(key=species,value=N,-c(time,grwr)) %>%
	mutate(species = recode(species, "Ev" = "Native", "Microstegium" = "Invader", "Microstegium litter" = "Invader litter"),
	time = time + 400)

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("BothInfection_PlantDiseaseModel_080519.pdf",width=3,height=3)
inf_low_tol_both_2 %>%
	filter(!(species %in% c("EvS", "EvA")) & time <= 600) %>%
	ggplot(aes(x = time, y = N, color = species)) +
	geom_line(size = 0.5)+
	scale_colour_manual(values=c("black", "gray", "#407879"), name = "") +
	theme_bw() +
	theme(axis.text=element_text(size=axisText,colour="black"),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position="none",
	plot.title=element_text(size=axisTitle)) +
	xlab("Years") +
	ylab("Population size")
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("BothInfection_NoAnn_PlantDiseaseModel_080519.pdf",width=3,height=3)
inf_low_tol_both_2 %>%
	filter(!(species %in% c("EvS", "EvA", "Invader")) & time <= 600) %>%
	ggplot(aes(x = time, y = N, color = species)) +
	geom_line(size = 1)+
	scale_colour_manual(values=c("gray", "#407879"), name = "") +
	theme_bw() +
	theme(axis.text=element_text(size=axisText,colour="black"),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position="none",
	plot.title=element_text(size=axisTitle)) +
	xlab("Years") +
	ylab("Population size")
dev.off()

filter(inf_low_tol_both, time == years)


## Sensitivity analysis function
sensFun=function(sensVal,paramNames,params,invader,N0.a,N0.p){

	# run model with default parameters
out_def=simFun(params=params,invader=invader,N0.a=N0.a,N0.s=0,N0.p=N0.p,L0=0,simtime=500)

	# save default output
	grwr_def=mean(filter(out_def, time <=400)$grwr)
	
	# adjust parameter values
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
  		out_star=simFun(params=p,invader=invader,N0.a=N0.a,N0.s=0,N0.p=N0.p,L0=0,simtime=500) 
  		
  		# grwr with perturbed values
  		grwr_star[i]=mean(filter(out_star, time <=400)$grwr)
  		
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
pinv=sensFun(sensVal=1.05,paramNames=paramNames,c(params, tolMin, tolMin),invader="P",N0.a=1,N0.p=0)
pinv

# Annual invasion
ainv=sensFun(sensVal=1.05,paramNames=paramNames,params=c(params, tolMin, tolMin),invader="A",N0.a=0,N0.p=1)
ainv

# Decrease by 5%

# Perennial invasion
pinv2=sensFun(sensVal=0.95,paramNames=paramNames,params=c(params, tolMin, tolMin),invader="P",N0.a=1,N0.p=0)
pinv2

# Annual invasion
ainv2=sensFun(sensVal=0.95,paramNames=paramNames,params=c(params, tolMin, tolMin),invader="A",N0.a=0,N0.p=1)
ainv2

# Combine dataframes
pinv$invader="native GRWR"
ainv$invader="invader GRWR"
pinv2$invader="native GRWR"
ainv2$invader="invader GRWR"

pinv$adjustment="increase 5%"
ainv$adjustment="increase 5%"
pinv2$adjustment="decrease 5%"
ainv2$adjustment="decrease 5%"

invdf=rbind(pinv,ainv,pinv2,ainv2)

# Absolute value of sensitivity
invdf$sensitivityAbs=abs(invdf$sensitivity)

# Plot
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("FullSensitivityAnalysis_PlantDiseaseModel_080519.pdf")
ggplot(invdf,aes(x=param,y=sensitivity))+geom_bar(stat="identity")+facet_grid(adjustment~invader)+geom_hline(yintercept=0,size=0.3)+theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("SensitivityAnalysis_PlantDiseaseModel_080519.pdf", width = 7, height = 6)
ggplot(filter(invdf,adjustment=="increase 5%"),aes(x=param,y=sensitivity))+
geom_bar(stat="identity",aes(fill=invader))+
facet_wrap(~invader)+
geom_hline(yintercept=0,size=0.3)+
scale_fill_manual(values=c("black", "#407879"))+
theme_bw() +
	theme(axis.text.y=element_text(size=(axisText-3),colour="black"),
	axis.text.x=element_text(size=(axisText-3),colour="black", angle = 45,hjust=1),
	axis.title=element_text(size=axisTitle,colour="black"),
	panel.background=element_blank(),
	panel.grid.major=element_blank(),
	panel.grid.minor=element_blank(),
	axis.line=element_line(colour="black"),
	panel.border=element_rect(color="black",fill=NA),
	legend.text=element_text(size=legendText),
	legend.title=element_text(size=legendTitle),
	legend.key=element_blank(),
	legend.position="none",
	plot.title=element_text(size=axisTitle),
	strip.text=element_text(size=axisTitle),
	strip.background=element_blank())+
	xlab("parameter")
dev.off()

setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
write.csv(invdf,"SensitivityAnalysis_PlantDiseaseModel_080519.csv",row.names=F)