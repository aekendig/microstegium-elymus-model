# clear all existing data
rm(list=ls())

# initial mass
init=6

# proportion remaining
c_rem_M=c(1,0.95,0.59,0.58,0.59)
c_rem_A=c(1,0.81,0.59,0.61,0.62)
e_rem_M=c(1,0.88,0.59,0.49,0.58)
e_rem_A=c(1,0.86,0.41,0.48,0.45)

# time
months=c(0,2,7,9,12)

# build dataframe
dat=data.frame(site=rep(rep(c("M","A"),each=5),2),invaded=rep(c(0,1),each=10),time=rep(months,4),mass.g=init*c(c_rem_M,c_rem_A,e_rem_M,e_rem_A),perc=c(c_rem_M,c_rem_A,e_rem_M,e_rem_A))
dat$invaded=as.factor(dat$invaded)

# plot
ggplot(dat,aes(x=time,y=mass.g))+geom_point(aes(colour=invaded,shape=site))+geom_line(aes(colour=invaded,linetype=site))

# overall mean decomposition after one year
mean(subset(dat,time==12)$mass.g)
mean(subset(dat,time==12)$perc)
