# Ofu Island ch1 CBASS - written by Courtney Klepac
# Physiological measurements - Fv/Fm from 0 to 21hrs##
setwd("~/Documents/Dissertation/Writing/manuscripts/Ch1/data/fvfm")
library(Hmisc)
library(MASS)
library(ggplot2)
library(lme4)
library(FSA)
library(qqplotr)
library(rcompanion)
require(multcomp)
require(nlme)
require(graphics)
library(stats)
library(dplyr)
library(qqplotr)
library(cowplot)
library(car)
library(psych)
library(PerformanceAnalytics)
library(tidyverse)

#read in file
pam<-read.csv("AS_PAM_ch1-2.csv", header=T)
#set factor levels
pam$origin=factor(pam$origin,levels(pam$origin)[c(1,3,2)])
pam$dest=factor(pam$dest,levels(pam$dest)[c(1,3,2)])
#pam$tank=factor(pam$tank,levels=unique(pam$tank))
pam$trt=factor(pam$trt,levels=(pam$trt)[c(2,1)])
pam$origin_dest=factor(pam$origin_dest,levels(pam$origin_dest)[c(1,4,2,5,3)])

#make coral colony discrete
pam$colony<-paste0(pam$colony,pam$origin)
as.factor(pam$colony)

summary(pam)
names(pam)

#creating interactions for post hoc comparisons
pam$grtrt<-interaction(pam$origin_dest,pam$trt)
pam$oritrt<-interaction(pam$origin,pam$trt)
pam$timetrt<-interaction(pam$time,pam$trt)
pam$grtime<-interaction(pam$origin_dest,pam$time)
pam$all<-interaction(pam$time,pam$origin_dest,pam$trt)

####Porites####
por<-pam[pam$species=='por',]
summary(por)
porheat<-por[por$trt=="heat",]
porcont<-por[por$trt=="cont",]
por$mqy<-(1-por$X21)

##testing normalized loss in fv/fm##
#looking for outliers
boxplot(Ylossnorm~origin_dest + time + trt, xlab="site", data=por)
#shapiro
aggregate(Ylossnorm ~ time + origin + dest +trt, data=por, FUN=function(x) shapiro.test(x)$p.value)
#bartlett HOV
for(i in c("16-Jan","16-Jul")){
	print(paste0("Timepoint",i))
	print(bartlett.test(Ylossnorm ~ origin_dest, data=por[por$time==i,]))
}
plotNormalHistogram(por$Ylossnorm)

resid<-lm(Ylossnorm~origin*dest*time, data=por)
plot(resid,which=1) #want random scatter; no apparent trendline
qqPlot(resid,data=por,na.action=na.exclude,envelope=0.95) #try a transform if data are non-normal


####repeated measures AOV -> to build in colony as a random effect####
###FIRST KEEP ORI AND DEST SEPARATE###
Ylossnorm.aov<-aov(Ylossnorm ~ time*origin*dest*trt + Error(colony), data=por)
summary(Ylossnorm.aov)
###Tukeys lme###
#origin 
model <- lme(Ylossnorm ~ origin, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(origin="Tukey")), test = adjusted(type = "bonferroni"))
#orign*trt
model <- lme(Ylossnorm ~ oritrt, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(oritrt="Tukey")), test = adjusted(type = "bonferroni"))
#trt
model <- lme(Ylossnorm ~ trt, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(trt="Tukey")), test = adjusted(type = "bonferroni"))

###COMBINING ORI AND DEST###
Ylossnorm.aov<-aov(Ylossnorm ~ time*origin_dest*trt + Error(colony), data=por)
summary(Ylossnorm.aov)

#time*ori_dest*trt
model <- lme(Ylossnorm ~ all, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(all="Tukey")), test = adjusted(type = "bonferroni"))
#origin_dest
model <- lme(Ylossnorm ~ origin_dest, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(origin_dest="Tukey")), test = adjusted(type = "bonferroni"))
#origin_dest*trt
model <- lme(Ylossnorm ~ grtrt, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(grtrt="Tukey")), test = adjusted(type = "bonferroni"))
#time
model <- lme(Ylossnorm ~ time, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(time="Tukey")), test = adjusted(type = "bonferroni"))


####Goniastrea####
gon<-pam[pam$species=='gon',]
summary(gon)
gonheat<-gon[gon$trt=="heat",]
goncont<-gon[gon$trt=="cont",]
gon$mqy<-(1-gon$X21)

##testing normalized loss in fv/fm##
#looking for outliers
boxplot(Ylossnorm~origin_dest + time + trt, xlab="site", data=gon)
#shapiro
aggregate(Ylossnorm ~ time + origin + dest +trt, data=gon, FUN=function(x) shapiro.test(x)$p.value)
#bartlett HOV
for(i in c("16-Jan","16-Jul")){
  print(paste0("Timepoint",i))
  print(bartlett.test(Ylossnorm ~ origin_dest, data=gon[gon$time==i,]))
}
plotNormalHistogram(gon$Ylossnorm)

resid<-lm(Ylossnorm~origin*dest*time, data=gon)
plot(resid,which=1) #want random scatter; no apparent trendline
qqPlot(resid,data=gon,na.action=na.exclude,envelope=0.95) #try a transform if data are non-normal


####repeated measures AOV -> colony as a random effect####
###FIRST KEEP ORI AND DEST SEPARATE###
Ylossnorm.aov<-aov(Ylossnorm ~ time*origin*dest*trt + Error(colony), data=gon)
summary(Ylossnorm.aov)
###Tukeys lme###
#time
model <- lme(Ylossnorm ~ time, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(time="Tukey")), test = adjusted(type = "bonferroni"))
#time*trt
model <- lme(Ylossnorm ~ timetrt, random = ~ 1 | colony, gon, na.action=na.exclude)
summary(glht(model, linfct=mcp(timetrt="Tukey")), test = adjusted(type = "bonferroni"))

##COMBINING ORIGIN AND DEST###
Ylossnorm.aov<-aov(Ylossnorm ~ time*origin_dest*trt + Error(colony), data=gon)
summary(Ylossnorm.aov)
#time 
model <- lme(Ylossnorm ~ time, random = ~ 1 | colony, gon, na.action=na.exclude)
summary(glht(model, linfct=mcp(time="Tukey")), test = adjusted(type = "bonferroni"))
#orign_dest*trt
model <- lme(Ylossnorm ~ grtrt, random = ~ 1 | colony, gon, na.action=na.exclude)
summary(glht(model, linfct=mcp(grtrt="Tukey")), test = adjusted(type = "bonferroni"))
#time*trt
model <- lme(Ylossnorm ~ timetrt, random = ~ 1 | colony, gon, na.action=na.exclude)
summary(glht(model, linfct=mcp(timetrt="Tukey")), test = adjusted(type = "bonferroni"))

##To generate a priori contrasts in an unbalanced random mixed model
miss1_mod<-lmer(Ylossnorm ~ all + (1 | colony), gon, na.action=na.exclude)
anova(miss1_mod)
summary(glht(miss1_mod,linfct=mcp(all="Tukey")), test = adjusted(type = "bonferroni"))

####Plots####
#Porites#
#Plot of 21hr values
#plot with dest as x and points are origin
sum=Summarize(X21~grtrt+time, data=por, digits=3)
sum$se = sum$sd / sqrt(sum$n)
sum$se = signif(sum$se, digits=3)
sum$Origin<-factor(c("HV","MV","LV","MV","LV","HV","MV","LV","MV","LV","HV","MV","LV","MV","LV","HV","MV","LV","MV","LV"))
#sum$Origin<-factor(c("HV","MV","LV","MV","LV","HV","MV","LV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","LV","HV","HV","HV","MV","LV","HV","HV","HV","MV","LV","HV","HV","HV","MV","LV"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$dest=factor(sum$dest,levels=unique(sum$dest))
sum$time=factor(sum$time,levels=unique(sum$time))
sum$trt=factor(c("cont","cont","cont","cont","cont","heat","heat","heat","heat","heat","cont","cont","cont","cont","cont","heat","heat","heat","heat","heat"))
sum$trt=factor(sum$trt,levels=unique(sum$trt))

#poster and publication figure
pd=position_dodge(.75)
plot.pam<-ggplot(sum,aes(x=dest,y=mean,fill=Origin,color=trt)) + 
  geom_point(aes(shape=Origin),size=3,position=pd) + 
  geom_errorbar(aes(ymin=mean - se,ymax=mean + se), width=.2,position=pd) + 
  facet_wrap(~time) +
  theme_bw() + 																		
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) + 
  scale_color_manual(values = c('gray','black'), name = "Treatment") +
  coord_cartesian(xlim=c(0.5,3.4), ylim=c(0,0.6), expand=F) + 
  ylab("Maximum Quantum Yield (Fv/Fm)\n after Acute Heat Stress") + 
  xlab("Transplant Site")
save_plot("20200214_por_fvfm21hrse.pdf", plot.pam,
          base_aspect_ratio = 1.3)

#boxplot to incorporate both heat and control 
pam <- ggplot(data=por, 
              aes(x=dest, y=mqy, label= time, fill=origin, color=trt)) +
  scale_fill_manual(values = c ("red", "gold", "blue"), name = "Reef site") +
  scale_color_manual(values = c('gray','black'), name = "Treatment") +
  stat_boxplot(geom ='errorbar', width = 0.5, lwd=0.5)+
  geom_boxplot(width=0.5, lwd=0.5) +
  expand_limits(y = 0)+
  facet_grid(~time, space = "free", scales = "free")+ #this can also be used but rows and columns have to be specified -> facet_grid(. ~ experiment)
  theme_bw() +
  xlab(label = "Transplant Site") + ylab(label = "Fv/Fm")
ggsave("./fvfm-trt_boxplot_20200209.pdf", width = 10, height = 6)

#Goniastrea#
summary(gon)
#plot with dest as x and points are origin
sum=Summarize(X21~grtrt+time, data=gon, digits=3)
sum$se = sum$sd / sqrt(sum$n)
sum$se = signif(sum$se, digits=3)
sum$Origin<-factor(c("HV","MV","LV","MV","HV","MV","LV","MV","HV","MV","LV","MV","LV","HV","MV","LV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","HV","HV","HV","MV","HV","HV","HV","MV","LV","HV","HV","HV","MV","LV"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$dest=factor(sum$dest,levels=unique(sum$dest))
sum$time=factor(sum$time,levels=unique(sum$time))
sum$trt=factor(c("cont","cont","cont","cont","heat","heat","heat","heat","cont","cont","cont","cont","cont","heat","heat","heat","heat","heat"))
sum$trt=factor(sum$trt,levels=unique(sum$trt))

plot.pamg<-ggplot(sum,aes(x=dest,y=mean,color=trt)) + 
  geom_point(aes(shape=Origin),size=4,position=pd) + 
  geom_errorbar(aes(ymin=mean - se,ymax=mean + se), width=.2,position=pd) + 
  facet_wrap(~time) +
  theme_bw() + 																		
  theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) + 
  scale_color_manual(values = c('gray','black'), name = "Treatment") +
  coord_cartesian(xlim=c(0.5,3.25), ylim=c(0,0.65), expand=F) + 
  ylab("Maximum Quantum Yield (Fv/Fm)\n after Acute Heat Stress") + 
  xlab("Transplant Site")
save_plot("20200214_gon_fvfm21hrse.pdf", plot.pamg,
          base_aspect_ratio = 1.3)
pdf(file="gonheat-X21-sd-panel-siteshape.pdf")
pd=position_dodge(.6)
ggplot(sum,aes(x=dest,y=mean,color=Origin)) + 
	geom_point(aes(shape=Origin),size=4,position=pd) + 
	geom_errorbar(aes(ymin=mean - sd,ymax=mean + sd), width=.2,position=pd) + 
	facet_wrap(~time) +
	theme_bw() + 																		
	theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) + 	scale_colour_manual(values=c("red","gold","blue")) + 						
	coord_cartesian(xlim=c(0.5,3.25), ylim=c(0,1.0), expand=F) + 
	ylab("Proportion of Maximum Quantum Yield (Fv/Fm) after Acute Heat Stress") + 
	xlab("Transplant Site")
dev.off()

#poster figure
pd=position_dodge(.75)
plot.pam<-ggplot(sum,aes(x=dest,y=mean,color=Origin)) + 
	geom_point(aes(shape=Origin),size=4,position=pd) + 
	geom_errorbar(aes(ymin=mean - sd,ymax=mean + sd), width=.2,position=pd) + 
	facet_wrap(~Time) +
	theme_bw() + 																		
	theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'), panel.grid.minor=element_blank()) + 	scale_colour_manual(values=c("red","gold","blue")) + 						
	coord_cartesian(xlim=c(0.5,3.25), ylim=c(0,1.0), expand=F) + 
	ylab("Proportion of Maximum Quantum Yield (Fv/Fm) after Acute Heat Stress") + 
	xlab("Transplant Site")
save_plot("20181220_gon_pam_poster.pdf", plot.pam,
          base_aspect_ratio = 1.3)
	
