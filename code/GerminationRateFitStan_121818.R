# clear all existing data
rm(list=ls())

# open libraries
library(plyr)
library(ggplot2)
library(Hmisc)

# plotting parameters
axisText=14
axisTitle=16
legendText=8
legendTitle=8
labelText=2
statText=2.5

# run "before-chapter.R"
source("/Users/AmyKendig/Google Drive/R Stuff/stan-book-master/before-chapter.R")

# set working directory
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/REU/Data")

# import data (all data combined up to post harvest)
dat=read.csv("MvEv_GerminationInfection_PostHarvest.csv",header=T)
head(dat)

#format date column
dat$Date=as.Date(as.character(dat$Date),"%Y%m%d") #ignore warning about time zone

#give records a pseudo date when they were counted very close to one another
ddply(dat,.(Date),summarise,nSamps=length(Date))

dat$Date[dat$Date=="2018-07-21"|dat$Date=="2018-07-22"|dat$Date=="2018-07-23"]="2018-07-24"
dat$Date[dat$Date=="2018-07-31"|dat$Date=="2018-08-01"]="2018-08-02"
dat$Date[dat$Date>="2018-08-22"]="2018-10-05"

ddply(dat,.(Date),summarise,nSamps=length(Date))

#days since planting
dat$Days=difftime(dat$Date,"2018-06-15",units="days")

#a numeric version of days
dat$Time=as.numeric(dat$Days)

#make a shade column
dat$Shade=with(dat,ifelse(Treatment=="Shade","yes","no"))

#make a column of litter levels (use treatment, but replace shade with medium)
dat$Litter=revalue(dat$Treatment,c("Shade"="Med"))
head(subset(dat,select=c(Treatment,Litter,Shade)),n=15)

#re-order litter levels
dat$Litter=factor(dat$Litter,levels=c("None","Low","Med","High"))

#quantitative litter levels
litterDat=data.frame(Litter=c("None","Low","Med","High"),Litter.g=c(0,0.91,1.82,3.64))
litterDat
dat2=merge(dat,litterDat,all=T)
head(dat2)

#re-order sp present levels
dat2$SpPresent=factor(dat2$SpPresent,c("Ev","Mv","Ev+Mv"))

#make a treatment ID
dat2$TrtID=with(dat2,paste(SpPresent,Litter,Shade,PotID,sep="."))

## Mv germination from litter

#Lili removed litter-derived Microstegium on 7/4 and 7/11. For 7/11, she included the number that were removed on 7/4. All counts after that were the observed numbers. Update new counts to include the plants that had been removed

#get the older data
mvTot=subset(subset(dat2,Litter%in%c("Low","Med","High")&SpPresent=="Ev"&Date=="2018-07-11"),select=c(TrtID,GermMv))
mvTot=plyr::rename(mvTot,c("GermMv"="GermMv2"))
mvTot

#select the newer dataset that needs to be updated
mvMiss=subset(subset(dat2,Litter%in%c("Low","Med","High")&SpPresent=="Ev"&(Date=="2018-07-24")),select=c(Date,TrtID,GermMv))
mvMiss

#merge older and newer data
newDat=merge(mvMiss,mvTot)
newDat
newDat$GermMv3=with(newDat,GermMv+GermMv2)

#use as older data
mvTot2=subset(newDat,select=c(TrtID,GermMv3))

#select newer data
mvMiss2=subset(subset(dat2,Litter%in%c("Low","Med","High")&SpPresent=="Ev"&(Date=="2018-08-02")),select=c(Date,TrtID,GermMv))
mvMiss2

#merge older and newer
newDat2=merge(mvMiss2,mvTot2)
newDat2
newDat2$GermMv4=with(newDat2,GermMv+GermMv3)

#use as older data
mvTot3=subset(newDat2,select=c(TrtID,GermMv4))

#select newer data
mvMiss3=subset(subset(dat2,Litter%in%c("Low","Med","High")&SpPresent=="Ev"&(Date=="2018-10-05")),select=c(Date,TrtID,GermMv))
mvMiss3 #all are NA

#merge older and newer
newDat3=merge(mvMiss3,mvTot3)
newDat3
newDat3$GermMv5=newDat3$GermMv4

#merge new data
head(newDat)
head(newDat2)
head(newDat3)

newDats=subset(newDat,select=-c(GermMv2))
newDats=plyr::rename(newDats,c("GermMv3"="NewGermMv"))
newDats
newDat2s=subset(newDat2,select=-c(GermMv3))
newDat2s=plyr::rename(newDat2s,c("GermMv4"="NewGermMv"))
newDat2s
newDat3s=subset(newDat3,select=-c(GermMv4))
newDat3s=plyr::rename(newDat3s,c("GermMv5"="NewGermMv"))
newDat3s

newDatTot=rbind(newDats,newDat2s,newDat3s)
newDatTot

#merge with main data
dat3=merge(dat2,newDatTot,all=T)
nrow(dat3)-nrow(dat2)
head(dat3)

#fill in NA's with germination data
dat3$NewGermMv[is.na(dat3$NewGermMv)]=dat3$GermMv[is.na(dat3$NewGermMv)]

#mean of mv germination in all treatments
corSum=ddply(dat3,.(Days,Litter,SpPresent),summarise,corGerm=round(mean(NewGermMv,na.rm=T))) # rounded means so that they're integers
head(corSum) 

#mean of mv germination in ev treatments
corSub=subset(subset(corSum,SpPresent=="Ev"),select=c(Days,Litter,corGerm))
corSub

#make the NA's from the none treatment into 0
corSub$corGerm[is.na(corSub$corGerm)]=0

#see why the value drops for Low
subset(dat3,SpPresent=="Ev"&Litter=="Low"&(Time==10|Time==12)) 
#Low 1 is at 3 on 6/25, then drops
#Low 2 is at 1 on 6/25, then drops
#These may be due to observational errors or plants dying. Either way, the same error or change could have occurred in other pots, so keep these numbers the same.

#merge data
dat4=merge(dat3,corSub,all=T)
head(dat4)


## Accidentally removed plants
#When Amy counted in July (7/24) and Lily counted in August (8/2) they accidentally removed some of the plants, which need to be added back

#subset notes
pulledNotes=unlist(lapply(dat4$Notes,grep,pattern="Pulled",ignore.case=T,value=T))
pulledNotes
pulledDat=subset(subset(dat4,(Date=="2018-07-24"|Date=="2018-08-02")&Notes%in%pulledNotes),select=c(Date,TrtID,Notes))
pulledDat

#add in the missing Ev
pulledDat$PulledEv=c(2,0,1,1,1,2,1,1,1,1,0,1,1) #assuming that when Lily didn't indicate whether it was Ev or Mv, it was Ev
pulledDat$PulledMv=c(0,0,0,0,0,1,0,0,0,0,1,0,0)

#separate by date
julPulledDat=subset(subset(pulledDat,Date=="2018-07-24"),select=-c(Date,Notes))
augPulledDat=subset(subset(pulledDat,Date=="2018-08-02"),select=-c(Date,Notes))

#update pulls from July
upDatAug=subset(dat4,Date>"2018-07-24")
upDatAug2=merge(upDatAug,julPulledDat,all=T)
head(upDatAug2)
upDatAug2$PulledEv[is.na(upDatAug2$PulledEv)]=0
upDatAug2$PulledMv[is.na(upDatAug2$PulledMv)]=0
upDatAug2$NewGermMv2=with(upDatAug2,NewGermMv+PulledMv)
upDatAug2$NewGermEv=with(upDatAug2,GermEv+PulledEv)

# update pulls from August
upDatSep=subset(subset(upDatAug2,select=-c(PulledEv,PulledMv)),Date>"2018-08-02")
upDatSep2=merge(upDatSep,augPulledDat,all=T)
head(upDatSep2)
upDatSep2$PulledEv[is.na(upDatSep2$PulledEv)]=0
upDatSep2$PulledMv[is.na(upDatSep2$PulledMv)]=0
upDatSep2$NewGermMv2=with(upDatSep2,NewGermMv2+PulledMv)
upDatSep2$NewGermEv=with(upDatSep2,NewGermEv+PulledEv)

#subset new data
upDat1=subset(subset(upDatAug2,Date<="2018-08-02"),select=-c(PulledEv,PulledMv))
head(upDat1)
upDat2=subset(upDatSep2,select=-c(PulledEv,PulledMv))
head(upDat2)

#merge update dat
upDat=merge(upDat1,upDat2,all=T)

#merge with main data
dat5=merge(dat4,upDat,all=T)
head(dat5)
dim(dat4)
dim(dat5)

#make the NA's into previous germ numbers
dat5$NewGermMv2[is.na(dat5$NewGermMv2)]=dat5$NewGermMv[is.na(dat5$NewGermMv2)]
dat5$NewGermEv[is.na(dat5$NewGermEv)]=dat5$GermEv[is.na(dat5$NewGermEv)]

#make corrections for litter-derived Mv 
dat5$CorGermMv=with(dat5,NewGermMv2-corGerm)
dat5$CorGermMv
#if values are negative, make these zero
dat5$CorGermMv[dat5$CorGermMv<0]=0

#check NA's
unique(subset(subset(dat5,is.na(CorGermMv)&SpPresent!="Ev"),select=c(Date,SpPresent,Notes)))
#either not to be used for germination or from the days when Ev only was counted (8/2)


## Remove pots not for germination
unique(dat5$Notes)
subset(dat5,Notes=="not to be used for germination data")
#didn't collect germination data on the dropped pot prior to 7/22
dat6=subset(dat5,TrtID!="Mv.Med.yes.3")
nrow(dat4)-nrow(dat5)


##Final formatting

#pot-scale infection prevalence (use observed total counts, note the Mv in Ev only pots on 7/11 include some non-observed values (pulled on 7/4))
dat6$propInfMv=with(dat6,InfectedMv/GermMv)
dat6$propInfEv=with(dat6,InfectedEv/GermEv)


## Summarize data

#summary by treatment
mvSum=ddply(subset(dat6,SpPresent!="Ev"),.(Date,Days,Time,Litter,Litter.g,Shade, SpPresent),summarise,meanGerm=mean((CorGermMv),na.rm=T),seGerm=sd(CorGermMv,na.rm=T)/sum(!is.na(CorGermMv)),meanGermUnCor=mean(NewGermMv2,na.rm=T),seGermUnCor=sd(NewGermMv2,na.rm=T)/sum(!is.na(NewGermMv2)),propInf=binconf(sum(InfectedMv,na.rm=T),sum(GermMv,na.rm=T))[1],lowerInf=binconf(sum(InfectedMv,na.rm=T),sum(GermMv,na.rm=T))[2],upperInf=binconf(sum(InfectedMv,na.rm=T),sum(GermMv,na.rm=T))[3])
mvSum
evSum=ddply(subset(dat6,SpPresent!="Mv"),.(Date,Days,Time,Litter,Litter.g, SpPresent),summarise,meanGerm=mean(NewGermEv,na.rm=T),seGerm=sd(NewGermEv,na.rm=T)/sum(!is.na(NewGermEv)),propInf=binconf(sum(InfectedEv,na.rm=T),sum(GermEv,na.rm=T))[1],lowerInf=binconf(sum(InfectedEv,na.rm=T),sum(GermEv,na.rm=T))[2],upperInf=binconf(sum(InfectedEv,na.rm=T),sum(GermEv,na.rm=T))[3],propInf2=mean(InfectedEv/GermEv,na.rm=T),seInf2=sd(InfectedEv/GermEv,na.rm=T)/sqrt(sum(!is.na(InfectedEv))))
evSum


## Decide which data to use
ggplot(subset(mvSum,!is.na(meanGerm)),aes(x=Days,y=meanGerm))+geom_point(aes(shape=Shade),size=2)+geom_line(aes(group=Shade),size=0.6)+geom_errorbar(width=0.1,aes(ymin=meanGerm-seGerm,ymax=meanGerm+seGerm,group=Shade))+facet_grid(SpPresent~Litter)+theme(strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA))+xlab("Days")+ylab("Mv germination number (50 planted)")
# are declines due to counting or germination in Ev pots?
ggplot(subset(mvSum,!is.na(meanGermUnCor)),aes(x=Days,y=meanGermUnCor))+geom_point(aes(shape=Shade),size=2)+geom_line(aes(group=Shade),size=0.6)+geom_errorbar(width=0.1,aes(ymin=meanGermUnCor-seGermUnCor,ymax=meanGermUnCor+seGermUnCor,group=Shade))+facet_grid(SpPresent~Litter)+theme(strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA))+xlab("Days")+ylab("Mv germination number (50 planted)")
# some yes, so use the last day
ggplot(evSum,aes(x=Days,y=meanGerm))+geom_point(size=2)+geom_line(size=0.6)+geom_errorbar(width=0.1,aes(ymin=meanGerm-seGerm,ymax=meanGerm+seGerm))+facet_grid(SpPresent~Litter)+theme(strip.background=element_blank(),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(color="black",fill=NA))+xlab("Days")+ylab("Ev germination number (50 planted)")
# Ev counts on last day are stems, don't use those, 6th day is almost always the max

unique(mvSum$Days)
unique(evSum$Days)

# subset data
mvDat=subset(dat6,Days>100&SpPresent=="Mv"&Shade=="no")
evDat=subset(dat6,round(Days,1)==38.8&SpPresent=="Ev")


## Plot litter effect
ggplot(mvDat,aes(x=Litter.g,y=CorGermMv))+geom_point()
ggplot(evDat,aes(x=Litter.g,y=NewGermEv))+geom_point()

# no value should be greater than 50 - reduce
mvDat$CorGermMv=ifelse(mvDat$CorGermMv>50,50,mvDat$CorGermMv)
evDat$NewGermEv=ifelse(evDat$NewGermEv>50,50,evDat$NewGermEv)

## Fit model

# change working directory
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R files/stan files")

# save data and run
mv_dat <- list(x=mvDat$Litter.g, y=mvDat$CorGermMv, n=rep(50,length(mvDat$CorGermMv)), J=length(mvDat$CorGermMv))
mv_germ_fit <- stan("germination.stan", data=mv_dat)
ev_dat <- list(x=evDat$Litter.g, y=evDat$NewGermEv, n=rep(50,length(evDat$NewGermEv)), J=length(evDat$NewGermEv))
ev_germ_fit <- stan("germination.stan", data=ev_dat)

# examine output
print(mv_germ_fit)
mv_simms=as.matrix(mv_germ_fit)
mv_g=median(mv_simms[,"g"])
mv_a=median(mv_simms[,"a"])

print(ev_germ_fit)
ev_simms=as.matrix(ev_germ_fit)
ev_g=median(ev_simms[,"g"])
ev_a=median(ev_simms[,"a"])

# plot
mvDat_pred=data.frame(Litter.g=seq(0,max(mvDat$Litter.g),length=20))
mvDat_pred$pred=mv_g/(1+mv_a*mvDat_pred$Litter.g)
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("GerminationRateFitStan_Mv_121818.pdf",width=5.8,height=5.8)
ggplot(mvDat,aes(x=Litter.g,y=CorGermMv/50))+geom_point(position=position_jitter(width=0.04))+geom_line(data=mvDat_pred,aes(x=Litter.g,y=pred))+ylim(0,1.01)+xlab("Litter (g)")+ylab("Proportion of Mv seeds germinated")+theme(axis.text=element_text(size=axisText,colour="black"),axis.title=element_text(size=axisTitle,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(linetype="solid",color="black",fill=NA,size=0.9),legend.text=element_text(size=legendText),legend.title=element_text(size=legendTitle),legend.key=element_blank(),legend.position="none")
dev.off()

evDat_pred=data.frame(Litter.g=seq(0,max(evDat$Litter.g),length=20))
evDat_pred$pred=ev_g/(1+ev_a*evDat_pred$Litter.g)
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("GerminationRateFitStan_Ev_121818.pdf",width=5.8,height=5.8)
ggplot(evDat,aes(x=Litter.g,y=NewGermEv/50))+geom_point(position=position_jitter(width=0.04))+geom_line(data=evDat_pred,aes(x=Litter.g,y=pred))+ylim(0,1.01)+xlab("Litter (g)")+ylab("Proportion of Ev seeds germinated")+theme(axis.text=element_text(size=axisText,colour="black"),axis.title=element_text(size=axisTitle,colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),panel.border=element_rect(linetype="solid",color="black",fill=NA,size=0.9),legend.text=element_text(size=legendText),legend.title=element_text(size=legendTitle),legend.key=element_blank(),legend.position="none")
dev.off()