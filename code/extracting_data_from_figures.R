#### Goal: extract data from figures

# reference: Emmanuel Jjunju, https://www.r-bloggers.com/digitizing-jpeg-graphs-in-r/


#### Set-up ####

# clear all existing data
rm(list=ls())

# libraries
library(jpeg) # to use jpeg images in R
library(zoo)

# digitize functions
ReadAndCal = function(fname){
  ReadImg(fname)
  calpoints <- locator(n=4,type='p',pch=4,col='blue',lwd=2)
  return(calpoints)
}

ReadImg = function(fname){
  img <- readJPEG(fname)
  op <- par(mar=c(0,0,0,0))
  on.exit(par(op))
  plot.new()
  rasterImage(img,0,0,1,1)
}

DigitData = function(col='red',type='p',...){
  type <- ifelse(type=='b','o',type)
  type <- ifelse(type%in%c('l','o','p'),type,'p')
  locator(type=type,col=col,...)
}

Calibrate = function(data,calpoints,x1,x2,y1,y2){
  x 		<- calpoints$x[c(1,2)]
  y 		<- calpoints$y[c(3,4)]
  
  cx <- lm(formula = c(x1,x2) ~ c(x))$coeff
  cy <- lm(formula = c(y1,y2) ~ c(y))$coeff
  
  data$x <- data$x*cx[2]+cx[1]
  data$y <- data$y*cy[2]+cy[1]
  
  return(as.data.frame(data))
}


#### Steps ####

# ReadAndCal opens the jpeg in a plotting window and lets you define points on the x and y axes. You must start by clicking on the left-most x-axis point, then the right-most axis point, followed by the lower y-axis point and finally the upper y-axis point. You don’t need to choose the end points of the axis, only two points on the axis that you know the x or y value for. As you click on each of the 4 points, the coordinates are saved in the object cal.

# DigitData returns you to the figure window, and now you can click on each of the data points you’re interested in retrieving values for. The function will place a dot (colored red in this case) over each point you click on, and the raw x,y coordinates of that point will be saved to the data.points list. When you’re finished clicking points, you need to hit stop/Finish or right-click to stop the data point collection.

# Calibrate converts those raw x,y coordinates into the same scale as the original graph. Feed the function your data.point list, the ReadAndCal list that contains your 4 control points from the first step, and then 4 numeric values that represent the 4 original points you clicked on the x and y axes. These values should be in the original scale of the figure (i.e. read the values off the graph’s tick marks).

#### Digitize figures ####

# Redwood et al. 2018 Fig. 2A
(cal_rw18_2a = ReadAndCal("./data/lit_figures/Redwood_2018_Fig2A.jpg"))
(data_rw18_2a = DigitData(col = 'red'))
df_rw18_2a = Calibrate(data_rw18_2a, cal_rw18_2a, 0, 23, 0, 100)
df_rw18_2a$x = round(df_rw18_2a$x)
df_rw18_2a$month = c("Dec 10", "Feb 11", "Apr 11", "Jul 11", "Sep 11", "Nov 11", "Jan 12", "Mar 12", "May 12", "Jul 12", "Sep 12", "Nov 12")
df_rw18_2a$surv = round(df_rw18_2a$y)

# Redwood et al. 2018 Fig. 2B
(cal_rw18_2b = ReadAndCal("./data/lit_figures/Redwood_2018_Fig2B.jpg"))
(data_rw18_2b = DigitData(col = 'red'))
df_rw18_2b = Calibrate(data_rw18_2b, cal_rw18_2b, 0, 23, 0, 1)
df_rw18_2b$x = round(df_rw18_2b$x)
df_rw18_2b$month = c("Dec 10", "Feb 11", "Apr 11", "Jul 11", "Sep 11", "Nov 11", "Jan 12", "Mar 12", "May 12", "Jul 12", "Sep 12", "Nov 12")
df_rw18_2b$germ = round(df_rw18_2b$y, 2)

# Foster and Gross 1998 Fig. 1B
(cal_fg98_1b = ReadAndCal("./data/lit_figures/Foster_1998_Fig1B.jpg"))
(data_fg98_1b = DigitData(col = 'red'))
df_fg98_1b = Calibrate(data_fg98_1b, cal_fg98_1b, 0, 374, 0, 1000)
df_fg98_1b$treatment = c("removed", "control", "added")
df_fg98_1b$litter.g = round(df_fg98_1b$y)

# Foster and Gross 1998 Fig. 3A
(cal_fg98_3a = ReadAndCal("./data/lit_figures/Foster_1998_Fig3A.jpg"))
(data_fg98_3a = DigitData(col = 'red'))
df_fg98_3a = Calibrate(data_fg98_3a, cal_fg98_3a, 0, 374, 0, 40)
df_fg98_3a$treatment = c("removed", "control", "added")
df_fg98_3a$germ = round(df_fg98_3a$y)/100

# Wilson et al. 2015 Fig. 3B
(cal_wi15_3b = ReadAndCal("./data/lit_figures/Wilson_2015_Fig3B.jpg"))
(data_wi15_3b = DigitData(col = 'red'))
df_wi15_3b = Calibrate(data_wi15_3b, cal_wi15_3b, 0, 40, 0, 20000)
df_wi15_3b$biomass.g = round(df_wi15_3b$x)
df_wi15_3b$seeds = round(df_wi15_3b$y)

# DeMeester and Richter 2010 Fig. 3
(cal_dr10_3 = ReadAndCal("./data/lit_figures/DeMeester_2010_Fig3.jpg"))
(data_dr10_3 = DigitData(col = 'red'))
df_dr10_3 = Calibrate(data_dr10_3, cal_dr10_3, 0, 400, 0, 120)
df_dr10_3$days = round(df_dr10_3$x)
df_dr10_3$mass.prop = round(df_dr10_3$y)/100
df_dr10_3$plot = c("M. vimineum", "weeded", "M. vimineum")

# Emery et al. 2013 Fig. 3A
(cal_em13_3A = ReadAndCal("./data/lit_figures/Emery_2013_Fig3.jpg"))
(data_em13_3A = DigitData(col = 'red'))
df_em13_3A = Calibrate(data_em13_3A, cal_em13_3A, 0, 1, 0, 700)
df_em13_3A$treatment = c(rep(c("reference", "spring09", "spring0910", "spring09herb"), each = 3), rep(c("spring10", "fall09", "fall09herb", "herb"), each = 2))
df_em13_3A$year = c(rep(c(2009, 2010, 2011), 4), rep(c(2010, 2011), 4))
df_em13_3A$seedlings = round(df_em13_3A$y)

# Emery et al. 2013 Fig. 3B
(cal_em13_3B = ReadAndCal("./data/lit_figures/Emery_2013_Fig3.jpg"))
(data_em13_3B = DigitData(col = 'red'))
df_em13_3B = Calibrate(data_em13_3B, cal_em13_3B, 0, 1, 0, 200)
df_em13_3B$treatment = c(rep(c("reference", "spring09", "spring0910", "spring09herb"), each = 2), "spring10", "fall09", "fall09herb", "herb")
df_em13_3B$year = c(rep(c(2009, 2010), 4), rep(2010, 4))
df_em13_3B$adults = round(df_em13_3B$y)

# combine Emery data
df_em13_3 <- merge(df_em13_3A[,c("treatment", "year", "seedlings")],
                   df_em13_3B[,c("treatment", "year", "adults")],
                   all = T)
df_em13_3$establishment = round(df_em13_3$adults / df_em13_3$seedlings, digits = 3)
plot(df_em13_3$seedlings, df_em13_3$establishment, type = "p")

# Warren et al. 2013 Fig. 2
(cal_wa13_2 = ReadAndCal("./data/lit_figures/Warren_2013_Fig2.jpg"))
(data_wa13_2 = DigitData(col = 'red'))
df_wa13_2 = Calibrate(data_wa13_2, cal_wa13_2, 0, 5, 0, 100)
df_wa13_2$litter_depth.cm = c(0, 1, 5)
df_wa13_2$seedling_survival = round(df_wa13_2$y/100, digits = 2)

# Reynolds et al. 2001 Fig. 1
(cal_ry01_1 = ReadAndCal("./data/lit_figures/Reynolds_2001_Fig1.jpg"))
(data_ry01_1 = DigitData(col = 'red'))
df_ry01_1 = Calibrate(data_ry01_1, cal_ry01_1, 0, 1, 0, 1)
df_ry01_1$litter = rep(c("bare", "light", "heavy"), each = 6)
df_ry01_1$species = rep(c("N. pulchra", "F. rubra", "C. nutkaensis", "D. holciformis", "F. arundinaceae", "H. lanatus"), 3)
df_ry01_1$germ.prop = round(df_ry01_1$y, digits = 2)

# Malmstrom et al. 2005 Fig. 2
(cal_ma05_2 = ReadAndCal("./data/lit_figures/Malmstrom_2005_Fig2.jpg"))
(data_ma05_2 = DigitData(col = 'red'))
df_ma05_2 = Calibrate(data_ma05_2, cal_ma05_2, 0, 1, 0, 250)
df_ma05_2$experiment = c(rep(1, 10), rep(2, 8))
df_ma05_2$species = c("EGM", "EGY", "EMC", "EMT", "HBY", "KML", "KMO", "KMY", "NPS", "NPY", "EEY", "EMC", "EMT", "NPC", "NPM", "NPS", "NPT", "NPY")
df_ma05_2$biomass.g = round(df_ma05_2$y, digits = 2)


#### Save values ####
write.csv(df_rw18_2a[,c("month", "surv")], "./data/Redwood_2018_Fig2A.csv", row.names = F)
write.csv(df_rw18_2b[,c("month", "germ")], "./data/Redwood_2018_Fig2B.csv", row.names = F)
write.csv(df_fg98_1b[,c("treatment", "litter.g")], "./data/Foster_1998_Fig1B.csv", row.names = F)
write.csv(df_fg98_3a[,c("treatment", "germ")], "./data/Foster_1998_Fig3A.csv", row.names = F)
write.csv(df_wi15_3b[,c("biomass.g", "seeds")], "./data/Wilson_2015_Fig3B.csv", row.names = F)
write.csv(df_dr10_3[,c("plot", "days", "mass.prop")], "./data/DeMeester_2010_Fig3.csv", row.names = F)
write.csv(df_em13_3, "./data/Emery_2013_Fig3.csv", row.names = F)
write.csv(df_wa13_2[,c("litter_depth.cm", "seedling_survival")], "./data/Warren_2013_Fig2.csv", row.names = F)
write.csv(df_ry01_1[,c("litter", "species", "germ.prop")], "./data/Reynolds_2001_Fig1.csv", row.names = F)
write.csv(df_ma05_2[,c("experiment", "species", "biomass.g")], "./data/Malmstrom_2005_Fig2.csv", row.names = F)


#### data calculations ####

# load packages
library(tidyverse)


#### convert litter depth to weight ####

# import Warren values
war <- read_csv("data/Warren_2013_Fig2.csv")

# values from Ash 1995
ash <- tibble(litter_weight.gm2 = c(544.76, 610.05, 471.30, 910.02, 544.57, 463.35, 524.94, 752.11),
              litter_depth.mm = c(62, 14, 64, 47, 58, 36, 60, 55)) %>%
  mutate(litter_depth.cm = litter_depth.mm * 0.1,
         cutting = rep(c("before", "after"), 4),
         cover = rep(rep(c("clear", "forest"), each = 2), 2),
         site = rep(c("BC", "RK"), each = 4),
         site_cover = paste(site, cover, sep = "_"))

# visualize
ggplot(ash, aes(litter_depth.cm, litter_weight.gm2)) +
  geom_point(aes(color = cutting, shape = site_cover)) +
  geom_smooth(method = "lm", formula = y ~ x)
# not a linear relationship

# regression
mod <- lm(litter_weight.gm2 ~ litter_depth.cm, data = ash)
summary(mod)

# average
mean(ash$litter_depth.cm) # 4.95
mean(ash$litter_weight.gm2) # 603

# calculate beta
war
# 0.08 = 0.96/(1 + B * 603)
# 1 + B * 603 = 0.96/0.08
# B * 603 = 0.96/0.08 - 1
(0.96/0.08 - 1)/603


#### perennial adult survival ####

# values from Lauenroth and Adler 2008
1 - 1/5
1 - 1/39


#### litter effects on perennials ####

# import data
rey <- read_csv("./data/Reynolds_2001_Fig1.csv")

# subset for native species
# select bare and high litter treatments
# make data wide
# caculate beta (equation in section above)
rey2 <- rey %>%
  filter(!(species %in% c("F. arundinaceae", "H. lanatus")) & litter != "light") %>%
  pivot_wider(names_from = litter,
              values_from = germ.prop) %>%
  mutate(beta = (bare/heavy - 1)/300)

range(rey2$beta)


#### perennial biomass range ####

# import data
malm <- read_csv("./data/Malmstrom_2005_Fig2.csv")
range(malm$biomass.g)
