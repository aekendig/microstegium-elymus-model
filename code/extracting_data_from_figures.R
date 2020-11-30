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
df_dr10_3 = Calibrate(data_dr10_3, cal_dr10_3, 0, 300, 0, 100)
df_dr10_3$days = round(df_dr10_3$x)
df_dr10_3$mass.prop = round(df_dr10_3$y)
df_dr10_3$plot = c("weeded", "M. vimineum")

# Emery et al. 2013 Fig. 3A
(cal_em13_3A = ReadAndCal("./data/lit_figures/Emery_2013_Fig3.jpg"))
(data_em13_3A = DigitData(col = 'red'))
df_em13_3A = Calibrate(data_em13_3A, cal_em13_3A, 0, 1, 0, 700)
df_em13_3A$treatment = c("reference", "reference")
df_em13_3A$year = c(2009, 2010)
df_em13_3A$seedlings = round(df_em13_3A$y)

# Emery et al. 2013 Fig. 3B
(cal_em13_3B = ReadAndCal("./data/lit_figures/Emery_2013_Fig3.jpg"))
(data_em13_3B = DigitData(col = 'red'))
df_em13_3B = Calibrate(data_em13_3B, cal_em13_3B, 0, 1, 0, 200)
df_em13_3B$treatment = c("reference", "reference")
df_em13_3B$year = c(2009, 2010)
df_em13_3B$adults = round(df_em13_3B$y)

# combine Emery data
df_em13_3 <- merge(df_em13_3A[,c("treatment", "year", "seedlings")],
                   df_em13_3B[,c("treatment", "year", "adults")],
                   all = T)
df_em13_3$establishment = df_em13_3$adults / df_em13_3$seedlings


#### Save values ####
write.csv(df_rw18_2a[,c("month", "surv")], "./data/Redwood_2018_Fig2A.csv", row.names = F)
write.csv(df_rw18_2b[,c("month", "germ")], "./data/Redwood_2018_Fig2B.csv", row.names = F)
write.csv(df_fg98_1b[,c("treatment", "litter.g")], "./data/Foster_1998_Fig1B.csv", row.names = F)
write.csv(df_fg98_3a[,c("treatment", "germ")], "./data/Foster_1998_Fig3A.csv", row.names = F)
write.csv(df_wi15_3b[,c("biomass.g", "seeds")], "./data/Wilson_2015_Fig3B.csv", row.names = F)
write.csv(df_dr10_3[,c("plot", "days", "mass.prop")], "./data/DeMeester_2010_Fig3.csv", row.names = F)
write.csv(df_em13_3, "./data/Emery_2013_Fig3.csv", row.names = F)
