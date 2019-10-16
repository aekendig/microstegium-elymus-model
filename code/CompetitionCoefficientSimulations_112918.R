# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(data.table)

# parameters
lambda=100 # fecundity
simtime=100 # simulation time
N=rep(NA,simtime) # population vector
N[1]=1

# intraspecific competition simulation
Nfun=function(lambda,simtime,N,alpha){
	for(i in 1:(simtime-1)){
		N[i+1]=N[i]*lambda/(1+alpha*N[i])
	}
	return(N)
}

# apply Nfun over a range of alpha values
alphas=c(0,0.1,0.5,1,5,10)
Nlist=list()
for(j in 1:length(alphas)){
	pop=Nfun(lambda,simtime,N,alphas[j])
	out=data.frame(time=1:simtime,pop=pop,pop2=c(pop[2:simtime],NA))
	out$cc=alphas[j]
	Nlist[[j]]=out
}

# dataframe
Ndf=data.frame(rbindlist(Nlist))
head(Ndf)
str(Ndf)
Ndf$cc=factor(Ndf$cc) #make cc a factor

# plot
ggplot(Ndf,aes(x=time,y=log(pop),colour=cc))+geom_line()

ggplot(subset(Ndf,cc!="0"),aes(x=time,y=pop,colour=cc))+geom_line()

ggplot(subset(Ndf,cc!="0"),aes(x=pop,y=pop2))+geom_line()+facet_wrap(~cc,scales="free")

# Final values
finN=subset(Ndf,time==simtime)
finN