#### set up ####

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

#### parameters ####

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
#alpha.sa=mean(c(0.724,9.574)) #5.149
#alpha.sa=9.574
alpha.sa=0.724
alpha.pa=alpha.sa/10 

#alpha.as=mean(c(0.054,0.357)) #0.2055
alpha.as=0.357
#alpha.ss=mean(c(0.002,0.006)) #0.004
alpha.ss=0.002
alpha.ps=0

#alpha.ap=mean(c(0.054,0.357)) #0.2055
alpha.ap=0.357
alpha.sp=alpha.ap 
alpha.pp=alpha.sp*10 

# convert seed production to biomass
c.a=30/6500

# decrease seed production due to infection
tolMin=0.19
tolMax=0.6

# parameter list
params=c(m.p,s.a,s.s,gamma.a,gamma.s,alpha.aL,alpha.sL,b,h.a,h.s,h.p,lambda.a,lambda.p,alpha.aa,alpha.as,alpha.ap,alpha.sa,alpha.ss,alpha.sp,alpha.pa,alpha.ps,alpha.pp,c.a)

# parameter names
paramNames=c("m.p","s.a","s.s","gamma.a","gamma.s","alpha.aL","alpha.sL","b","h.a","h.s","h.p","lambda.a","lambda.p","alpha.aa","alpha.as","alpha.ap","alpha.sa","alpha.ss","alpha.sp","alpha.pa","alpha.ps","alpha.pp","c.a","tol.a","tol.p")

# simulation time
years = 500


#### population model ####

simFun=function(params,N0.a,N0.s,N0.p,L0,Ni.a,simtime){
  
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
    
    # introduce annual
    N.a[t]=ifelse(N0.a==0 & t==100, Ni.a, N.a[t])
    
    # calulate parameters to introduce disease
    lam.a=ifelse(t<200, lambda.a, tol.a*lambda.a)
    lam.p=ifelse(t<200, lambda.p, tol.p*lambda.p)
    lam.s=ifelse(t<200, lambda.s, tol.p*lambda.s)
    
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
  }
  
  # save data
  dfN=data.frame(time=rep(1:simtime,4),N=c(N.s,N.p,N.a,L),species=rep(c("Elymus seedling","Elymus adult","Microstegium","Microstegium litter"),each=simtime))
  
  # return
  return(dfN)
}


#### no disease over time ####
simFun(params=c(params, 1, 1),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=years) %>%
  ggplot(aes(x = time, y = log(N+1), color = species)) +
  geom_line()

#### disease in annual over time ####
simFun(params=c(params, tolMin, 1),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=years) %>%
  ggplot(aes(x = time, y = log(N+1), color = species)) +
  geom_line()

#### disease in both over time ####
simFun(params=c(params, tolMin, tolMin),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=years) %>%
  ggplot(aes(x = time, y = log(N+1), color = species)) +
  geom_line()


#### final density over parameter ####

# define function
parmFun <- function(parm_ID, parm_min, parm_max, parm_length){
  
  # parameter values
  parm_vals = seq(parm_min, parm_max, length.out = parm_length)
  
  # store final values
  parm_dat = tibble(species = NA, N = NA, parm = NA)
  
  for(i in 1:parm_length) {
    
    # set parameter
    params_temp <- c(params, tolMin, tolMin)
    params_temp[parm_ID] <- parm_vals[i]
    
    # run simulation
    sim <- simFun(params=params_temp,N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=years)
    
    # summarise last 50 years
    sim_sum <- sim %>%
      filter(time >= (years - 50)) %>%
      group_by(species) %>%
      summarise(N = mean(N)) %>%
      ungroup()
    
    # add parameter values
    sim_sum$parm = parm_vals[i]
    
    # save last row
    parm_dat <- rbind(parm_dat, sim_sum)
  }
  
  # remove first 
  parm_dat <- filter(parm_dat, !is.na(parm))
  
  # plot
  print(ggplot(parm_dat, aes(x = parm, y = log(N+1), color = species)) +
    geom_line())
}

# lambda.a
parmFun(12, 0, 6000, 20)

# alpha.aa
parmFun(14, 0.001, 9, 20)

# alpha.aa smaller
parmFun(14, 0, 0.005, 20)
parmFun(14, 0, 0.01, 100)
parmFun(14, 0, 0.1, 100)
parmFun(14, 0, 0.5, 100)

# alpha.sa
parmFun(17, 0.001, 10, 20) # no change


#### different alpha.aa over time ####

# try value based on parmFun
parms_temp <- params
parms_temp[14] <- 0.15

# without disease
simFun(params=c(parms_temp, 1, 1),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=800) %>%
  ggplot(aes(x = time, y = log(N+1), color = species)) +
  geom_line()

# with disease
simFun(params=c(parms_temp, tolMin, tolMin),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=800) %>%
  ggplot(aes(x = time, y = log(N+1), color = species)) +
  geom_line()
# not sure if it's justified to multiply alpha.aa by 10


#### presentation figure ####

# healthy population
pres_dat_h <- simFun(params=c(parms_temp, 1, 1),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=800) %>%
  spread(.,species,N) %>% 
  rename(.,EvA="Elymus adult",EvS="Elymus seedling")%>% 
  mutate(Ev=EvA+EvS)%>% 
  gather(key=species,value=N,-time) %>%
  mutate(species = recode(species, "Ev" = "native", "Microstegium" = "invader"),
         scenario = "no disease")

# infected population
pres_dat_i <- simFun(params=c(parms_temp, tolMin, tolMin),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=800) %>%
  spread(.,species,N) %>% 
  rename(.,EvA="Elymus adult",EvS="Elymus seedling")%>% 
  mutate(Ev=EvA+EvS)%>% 
  gather(key=species,value=N,-time) %>%
  mutate(species = recode(species, "Ev" = "native", "Microstegium" = "invader"),
         scenario = "disease (native susceptible)")

# perennial tolerant
pres_dat_t <- simFun(params=c(parms_temp, tolMin, 1),N0.a=0,N0.s=1,N0.p=1,L0=0,Ni.a=1,simtime=800) %>%
  spread(.,species,N) %>% 
  rename(.,EvA="Elymus adult",EvS="Elymus seedling")%>% 
  mutate(Ev=EvA+EvS)%>% 
  gather(key=species,value=N,-time) %>%
  mutate(species = recode(species, "Ev" = "native", "Microstegium" = "invader"),
         scenario = "disease (native resistant)")


# combine data
pres_dat_nd <- pres_dat_h %>%
  filter(species %in% c("native", "invader")) %>%
  mutate(species = fct_relevel(species, "native"),
         scenario = fct_relevel(scenario, "no disease"))

pres_dat_tol <- pres_dat_t %>%
  filter(time >= 200) %>%
  full_join(pres_dat_h) %>%
  filter(species %in% c("native", "invader")) %>%
  mutate(species = fct_relevel(species, "native"),
         scenario = fct_relevel(scenario, "no disease"))

pres_dat <- pres_dat_i %>%
  full_join(pres_dat_t) %>%
  filter(time >= 200) %>%
  full_join(pres_dat_h) %>%
  filter(species %in% c("native", "invader")) %>%
  mutate(species = fct_relevel(species, "native"),
         scenario = fct_relevel(scenario, "no disease"))

# plot with no disease
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("Simulation_NoDis_PlantDiseaseModel_LitOnly_080819.pdf",width=9,height=5)
ggplot(pres_dat_nd, aes(x = time, y = log(N+1), color = species, linetype = scenario)) +
  geom_line(size = 1.2)+
  scale_colour_manual(values=c("#407879", "black"), name = "Species") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted"), name = "Scenario") +
  theme_bw() +
  theme(axis.text=element_text(size=axisText,colour="black"),
        axis.title=element_text(size=axisTitle,colour="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        panel.border=element_rect(color="black",fill=NA),
        legend.text=element_text(size=legendText),
        legend.title=element_text(size=axisText),
        legend.key=element_blank(),
        legend.position=c(0.7, 0.15),
        plot.title=element_text(size=axisTitle),
        legend.box="horizontal",
        legend.key.width = unit(1,"cm")) +
  xlab("Years") +
  ylab("ln(Population size)")
dev.off()

# plot with tolerant native
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("Simulation_TolNat_PlantDiseaseModel_LitOnly_080819.pdf",width=9,height=5)
ggplot(pres_dat_tol, aes(x = time, y = log(N+1), color = species, linetype = scenario)) +
  geom_line(size = 1.2)+
  scale_colour_manual(values=c("#407879", "black"), name = "Species") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted"), name = "Scenario") +
  theme_bw() +
  theme(axis.text=element_text(size=axisText,colour="black"),
        axis.title=element_text(size=axisTitle,colour="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        panel.border=element_rect(color="black",fill=NA),
        legend.text=element_text(size=legendText),
        legend.title=element_text(size=axisText),
        legend.key=element_blank(),
        legend.position=c(0.7, 0.15),
        plot.title=element_text(size=axisTitle),
        legend.box="horizontal",
        legend.key.width = unit(1,"cm")) +
  xlab("Years") +
  ylab("ln(Population size)")
dev.off()

# plot with all three
setwd("/Users/AmyKendig/Google Drive/Microstegium Bipolaris/R Output")
pdf("Simulation_All_PlantDiseaseModel_LitOnly_080819.pdf",width=9,height=5)
ggplot(pres_dat, aes(x = time, y = log(N+1), color = species, linetype = scenario)) +
  geom_line(size = 1.2)+
  scale_colour_manual(values=c("#407879", "black"), name = "Species") +
  scale_linetype_manual(values=c("solid", "dashed", "dotted"), name = "Scenario") +
  theme_bw() +
  theme(axis.text=element_text(size=axisText,colour="black"),
        axis.title=element_text(size=axisTitle,colour="black"),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        panel.border=element_rect(color="black",fill=NA),
        legend.text=element_text(size=legendText),
        legend.title=element_text(size=axisText),
        legend.key=element_blank(),
        legend.position=c(0.7, 0.15),
        plot.title=element_text(size=axisTitle),
        legend.box="horizontal",
        legend.key.width = unit(1,"cm")) +
  xlab("Years") +
  ylab("ln(Population size)")
dev.off()