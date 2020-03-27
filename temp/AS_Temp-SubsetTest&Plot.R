setwd("~/Documents/Dissertation/Writing/manuscripts/Ch1/data/temp/") 
getwd()
library(Hmisc)
library(MASS)
library(ggplot2)
library(lme4)
library(multcomp)
library(lubridate)
library(data.table)
require(lsmeans)
library(multcompView)
library(car)
library(dplyr)
library(tidyr)
library(Hmisc)
library(FSA)

#Ofu pool temperatures 2015-2016
#read in hobotemp file
temp<-read.delim("AStemp2016-17_avg2.txt",header=T)
temp<-read.delim(pipe("pbpaste"))
head(temp)
#strip Date and time
temp$DateTime<-strptime(temp$DateTime, format="%m/%d/%y %H:%M")
#make month and year-month column to do summarize stats
temp$year=year(temp$DateTime)
temp$month=month(temp$DateTime)
temp$day=day(temp$DateTime)
temp<-unite(temp, 'date', c('month','year'),remove=F)
temp<-unite(temp, 'd', c('day','month','year'),remove=F)

#temp=mutate(temp,date=factor(date,levels=unique(date)))
str(temp)
tail(temp)

#summarize max, min, mean, etc for each day
temp$d=factor(temp$d,levels=unique(temp$d))
HV<-Summarize(HV~d,data=temp,digits=3)
MV<-Summarize(MV~d,data=temp,digits=3)
LV<-Summarize(LV~d,data=temp,digits=3)
write.table(HV,"ASdailytempsumm2-HV_2016-17.txt",sep = "\t")
write.table(MV,"ASdailytempsumm2-MV_2016-17.txt",sep = "\t")
write.table(LV,"ASdailytempsumm2-LV_2016-17.txt",sep = "\t")

#summarize max, min, mean, etc for each month
HVmo<-Summarize(HV~date,data=temp,digits=3)
MVmo<-Summarize(MV~date,data=temp,digits=3)
LVmo<-Summarize(LV~date,data=temp,digits=3)

#black bgnd plot
pdf("2015-16_AStempplot-blackbg.pdf",11,7)
par(bg="black",fg="white",col.lab="white")
plot(temp$DateTime, temp$HV, type="l", col="red", lwd=1.5, ylab="Temperature (°C)", xlab="2015-2016", main="Ofu HV Pool Temperature",col.main="white",col.axis="white",xaxt='n',ylim=c(25.5,35))
points(temp$DateTime, temp$MV, type="l", col="gold", lwd=1.5)
points(temp$DateTime, temp$LV, type="l", col="blue", lwd=1.5)
abline(h=30.2,lty=2,col="white")
#text(0,30.5,"Bleaching Threshold 30.1 °C",col="black")
#to specify the x-axis by months within the DateTime, and remove tick marks
axis.POSIXct(side=1,at=temp$DateTime,format='%b-%y',col.axis="white",lwd=0,lwd.tick=0)
legend("topright",c("HV","MV","LV"),lty=1,lwd=2.5,col=c("red","gold","blue"),bty="n")
dev.off()

#white background
pdf("2015-16_AStempplot.pdf",11,7)
plot(temp$DateTime, temp$HV, type="l", col="red", lwd=1.5, ylab="Temperature (°C)", xlab="2015-2016", main="Ofu HV Pool Temperature",xaxt='n',ylim=c(25.5,35))
points(temp$DateTime, temp$MV, type="l", col="gold", lwd=1.5)
points(temp$DateTime, temp$LV, type="l", col="blue", lwd=1.5)
abline(h=30.2,lty=2,col="black")
#text(0,30.5,"Bleaching Threshold 30.1 °C",col="black")
#to specify the x-axis by months within the DateTime, and remove tick marks
axis.POSIXct(side=1,at=temp$DateTime,format='%b-%y',col.axis="black",lwd=0,lwd.tick=0)
legend("topright",c("HV","MV","LV"),lty=1,lwd=2.5,col=c("red","gold","blue"),bty="n")
dev.off()


#copied summary stats of each site into one data frame and saved as tempsumday.txt in Excel. Imported here.
temp<-read.delim("~/Documents/Dissertation/Writing/manuscripts/Ch1/data/temp/AS_2015-16_tempsumday.txt",header=T)
str(temp)
summary(temp)
temp$site=factor(temp$site,levels(temp$site)[c(1,3,2)])

#anova of max temp
model=lm(max~site*season,data=temp)
anova(model)
summary(model)
##
lsmeans(model,pairwise~site*season,adjust='tukey')
lsmeans(model,pairwise~season,adjust='tukey')

#anova of mean temp
model=lm(mean~site*season,data=temp)
anova(model)
summary(model)
##
lsmeans(model,pairwise~season,adjust='tukey')

#anova of min temp
model=lm(min~site*season,data=temp)
anova(model)
summary(model)
##
lsmeans(model,pairwise~site,adjust='tukey')
lsmeans(model,pairwise~season,adjust='tukey')

#plotting monthly summary of mean, max, min 
summary(temp)
temp$d<-as.Date(strptime(temp$d , "%m/%d/%y"))
temp$month <- as.factor(format(temp$d,'%b'))
class(temp$month)
levels(temp$month)
temp$month=factor(temp$month,levels(temp$month)[c(6,7,3,4,5,1,2,10,9,8,11,12)])

#boxplot of mean
pdf("2015-16_OfuDailyboxplot.pdf",10,8)

p1<-ggplot(temp,aes(x=month, y=mean, fill=site)) + geom_boxplot() +         
xlab("Month") + ylab("Temperature (°C)") + ggtitle("A. Mean Daily Temperatures")+ scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw() 

p2<-ggplot(temp,aes(x=month, y=max, fill=site)) + geom_boxplot() +         
xlab("Month") + ylab("Temperature (°C)") + ggtitle("B. Max Daily Temperatures") + scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw()  
 
p3<-ggplot(temp,aes(x=month, y=range, fill=site)) + geom_boxplot() +         
xlab("Month") + ylab("Temperature (°C)") + ggtitle("C. Daily Temperature Range")+ scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw() #+ scale_color_manual(name="Site",values=c("black","black","black")) + theme(legend.position="bottom")
   
source("http://peterhaschke.com/Code/multiplot.R")       
multiplot(p1, p2, p3, cols=1)       
dev.off()   

#seasonal summary of mean, max, min 
#boxplot of mean
library(cowplot)
pdf("20190304_2015-16_OfuSeasonboxplot-outliersmean.pdf",10,8)

p3<-ggplot(temp,aes(x=season, y=mean, fill=site)) + geom_boxplot() + stat_summary(fun.y=mean, geom="point", shape=18, size=3, position = position_dodge(0.75)) +       
 ylab("Mean Daily Temperature (°C)") + ggtitle("C") + scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw(base_size=14) +theme(axis.title.x=element_blank())

p1<-ggplot(temp,aes(x=season, y=max, fill=site)) + geom_boxplot() +  stat_summary(fun.y=mean, geom="point", shape=18, size=3, position = position_dodge(0.75)) +              
ylab("Max Daily Temperature (°C)") + ggtitle("A") + scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw(base_size=14) +theme(axis.title.x=element_blank())  
 
p2<-ggplot(temp,aes(x=season, y=range, fill=site)) + geom_boxplot() +  stat_summary(fun.y=mean, geom="point", shape=18, size=3, position = position_dodge(0.75)) +              
ylab("Daily Temperature Range (°C)") + ggtitle("B")+ scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw(base_size=14) +theme(axis.title.x=element_blank()) 
plot_grid(p1, p2, p3, nrow=1)    
#save_plot("20190221_2015-16_OfuSeasonboxplot-outliers-6.pdf", gg_all, base_height = 6)
dev.off()   
       
    
#summary of mean, max, min by site
#boxplot of mean
pdf("2015-16_OfuTempSummboxplot-outliers.pdf",10,8)

p3<-ggplot(temp,aes(x=site, y=mean, fill=site)) + geom_boxplot() +         
xlab("Site") + ylab("Temperature (°C)") + ggtitle("C. Mean Daily Temperatures")+ scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw() 

p1<-ggplot(temp,aes(x=site, y=max, fill=site)) + geom_boxplot() +         
xlab("Site") + ylab("Temperature (°C)") + ggtitle("A. Max Daily Temperatures") + scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw()  
 
p2<-ggplot(temp,aes(x=site, y=range, fill=site)) + geom_boxplot() +         
xlab("Site") + ylab("Temperature (°C)") + ggtitle("B. Daily Temperature Range")+ scale_fill_manual(guide=FALSE,values=c("red","gold","blue")) + theme_bw() #+ scale_color_manual(name="Site",values=c("black","black","black")) + theme(legend.position="bottom")
    
source("http://peterhaschke.com/Code/multiplot.R")       
multiplot(p1, p2, p3, cols=3)        
dev.off()       

#######################
#ANOVA of Ofu historical temperatures: DHW,SST,SSTA 
noaa<-read.csv("~/Documents/Dissertation/Writing/manuscripts/Ch1/data/temp/NOAACRW5km_2000-2017_OfuIsland.csv")
noaa=mutate(noaa,year=factor(year,levels=unique(year)))
noaa=mutate(noaa,month=factor(month,levels=unique(month)))
noaa=mutate(noaa,day=factor(day,levels=unique(day)))
noaa=mutate(noaa,month=month.abb[month])
str(noaa)
library(tidyr)
noaa<-unite(noaa, date, c(year,month),remove=F)
noaa=mutate(noaa,date=factor(date,levels=unique(date)))
summary(noaa)
#summary of each variable, separated by each month
Summarize(DHW~date,data=noaa,digits=3)
Summarize(SST~date,data=noaa,digits=3)
Summarize(SSTA~date,data=noaa,digits=3)
bartlett.test(DHW~date,data=noaa) #if pval >= 0.05, variance is equal

#anova of DHW
model=lm(DHW~year,data=noaa)
anova(model)
summary(model)

lsmeans(model,pairwise~year,adjust='tukey')
tukey<-cld(leastsquare,alpha=0.05,Letters=letters,adjust='tukey')
plot(tukey,las=1,col="black") 

#anova of SST
model=lm(SST~year,data=noaa)
anova(model)
lsmeans(model,pairwise~year,adjust='tukey')

#anova of SSTA
model=lm(SSTA~year,data=noaa)
anova(model)
lsmeans(model,pairwise~year,adjust='tukey')
 

########################
#line plot of temps from 2010-2012 & 2015-2017
stan<-noaa[noaa$year =='2010' | noaa$year == '2011' | noaa$year == '2012',]
odu<-noaa[noaa$year =='2015' | noaa$year == '2016'| noaa$year == '2017',]
write.csv(stan,"NOAA2010-12.csv")
write.csv(odu,"NOAA2015-17.csv")

dat<-read.csv("~/Documents/Dissertation/Writing/manuscripts/Ch1/data/temp/NOAAstanVme.csv")
summary(dat)
dat=mutate(dat,date=factor(date,levels=unique(date)))
#dat$date <- strptime(dat$date,format="%b-%d")
names(dat)
head(dat)
#Set the Theshold Above Which to Count Temperatures
Threshold=30.2

#plot of 2015-2017 vs. 2010-2012 SST & DHW
pdf("NOAA-stanVSme-SST&DHW.pdf",11,7)
par(mar=c(5,5,3,5))
#plot SST for each year
plot(dat$date, dat$SST10, type='n', cex.lab=1.5,cex.axis=1.5, ylab="Temperature (°C)", xlab="Month",main="Ofu Pool Temperature & DHW",ylim=c(24,32))
points(dat$date, dat$SST10, type="l", col="red", lwd=3)
points(dat$date, dat$SST11, type="l", col="orange", lwd=3)
points(dat$date, dat$SST12, type="l", col="green3", lwd=3)
points(dat$date, dat$SST15, type="l", col="blue", lwd=3)
points(dat$date, dat$SST16, type="l", col="steelblue1", lwd=3)
points(dat$date, dat$SST17, type="l", col="mediumpurple", lwd=3)
abline(h=30.2)
#axis(1,at=dat$date,tick=F,labels=format(as.Date(dat$date), "%b"))
#plot DHW for each year under SST, extend y axis to leave room for DHW
par(new=T)
plot(dat$date,dat[,7],type="l",lty=5, col="blue", lwd=3, xaxt='n',xlab='',yaxt='n',ylab='',ylim=c(0,20))
points(dat$date,dat[,15],type="l",lty=5, col="red", lwd=3)
points(dat$date,dat[,19],type="l",lty=5, col="orange", lwd=3)
points(dat$date,dat[,23],type="l",lty=5, col="green3", lwd=3)
points(dat$date,dat[,11],type="l",lty=5, col="steelblue1", lwd=3)
points(dat$date,dat[,27],type="l",lty=5, col="mediumpurple", lwd=3)
axis(side=4, cex.axis=1.5)
mtext("DHW (°C week)",cex=1.5,side=4,line=3)
legend(c(as.POSIXct("2018-12-26 EST"),5),cex=1.5,c("2010","2011","2012","2015","2016","2017"),lty=1,lwd=4,col=c("red","orange","green3","blue","steelblue1","mediumpurple"),bty="n")
dev.off()


#median box plot of historical DHW
pdf("2000-17_OfuDHW-medianboxplot.pdf",8,8)
ggplot(noaa, aes(x = year, y = DHW, group=year))  +
       geom_boxplot(colour = "black", width = 0.5)  +         
       labs(x = "Year", 
            y = "DHW ")  +
       ## ggtitle("Main title") + 
       theme_classic()  + theme(panel.grid.major.y = element_line( size=.1, color="gray"))
dev.off()


##frequency plot of daily temp ranges
pdf("OfuTempRange.pdf",5,5)
ggplot(temp)+geom_freqpoly(aes(x=range,color=site,linetype=site),size=1)+theme_classic()+xlab("Daily Temperature Range (°C)")+ylab("Frequency of Occurence (%)")+scale_color_manual(guide=FALSE,values=c("red","gold","blue"))
dev.off()
 
#box plot of median
plot(noaa$DHW,data=noaa)
library(ggplot2)
letters=c("a","bcd","e","f","ab","cde","ab","de","a","abc","de","a","ab","abcd","cde","g","f","h")
pdf("2000-17_OfuDHW-medianboxplot.pdf",8,8)
ggplot(noaa, aes(x = year, y = DHW, group=year))  +
       geom_boxplot(colour = "black", width = 0.5)  +         
       labs(x = "Year", 
            y = "DHW ")  +
       ## ggtitle("Main title") + 
       theme_classic()  + theme(panel.grid.major.y = element_line( size=.1, color="gray"))
dev.off()
             
#boxplot with the mean instead median
pdf("2000-17_OfuDHW.pdf",8,8)
ggplot(noaa, aes(year)) +
geom_boxplot(aes(ymin = min, lower = Q1, middle=mean, upper=Q3,ymax=max),stat="identity") +
	xlab("Year") + ylab("Degree Heating Weeks") +
	theme_classic() + #title("Mean Degree Heating Weeks")
	theme( # remove the vertical grid lines
           panel.grid.major.x = element_blank() ,
           # explicitly set the horizontal lines (or they will disappear too)
           panel.grid.major.y = element_line( size=.1, color="gray")) 
dev.off()
