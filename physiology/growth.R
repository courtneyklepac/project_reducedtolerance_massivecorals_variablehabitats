# Ofu Island ch1 CBASS - written by Courtney Klepac
# Physiological measurements - weekly growth rate
setwd("~/Documents/Dissertation/Writing/manuscript/Ch1/data/")
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
growth<-read.csv("AS_growth-repfix.csv", header=T)
#set factor levels
growth$colony<-paste0(growth$colony,growth$origin)
growth$origin=factor(growth$origin,levels=unique(growth$origin))
growth$dest=factor(growth$dest,levels=unique(growth$dest))
growth$origin_dest=factor(growth$origin_dest,levels(growth$origin_dest)[c(1,4,2,5,3)])
growth$origin_dest=factor(growth$origin_dest,levels=unique(growth$origin_dest))
growth$time=factor(growth$time,levels(growth$time)[c(2,1)])
#creating interactions for post hoc comparisons
growth$grtime<-interaction(growth$time,growth$origin_dest) 
growth$oritime<-interaction(growth$origin,growth$time)
growth$desttime<-interaction(growth$dest,growth$time)

summary(growth)
str(growth)

####Porites####
#making two data frames for ea species
por<-growth[growth$species=='por',]
#porwin<-growth[growth$time=='16-Jul',]
#porsum<-growth[growth$time=='16-Jan',]
head(por)
summary(por)

#looking for outliers
boxplot(grate ~ origin_dest + time, data=por, las=2)
#shapiro
aggregate(grate ~ time + origin + dest, data=por, FUN=function(x) shapiro.test(x)$p.value)
#bartlett HOV
for(i in c("16-Jan","16-Jul")){
	print(paste0("Timepoint",i))
	print(bartlett.test(grate ~ origin_dest, data=por[por$time==i,]))
}
plotNormalHistogram(por$grate)

resid<-lm(grate~origin_dest*time, data=por)
plot(residuals(resid)) #want random scatter; no apparent trendline
summary(resid)
qqPlot(resid,data=por) 


##UNTRANSFORMED W/ ORI AND DEST SEPARATE
grate.aov<-aov(grate~ time*origin*dest + Error(colony), data=por)
summary(grate.aov)
model <- lme(grate ~ origin, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(origin="Tukey")), test = adjusted(type = "bonferroni"))
#time
model <- lme(grate ~ time, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(time="Tukey")), test = adjusted(type = "bonferroni"))
#origin*time
model <- lme(grate ~ oritime, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(oritime="Tukey")), test = adjusted(type = "bonferroni"))
#dest
model <- lme(grate ~ dest, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(dest="Tukey")), test = adjusted(type = "bonferroni"))
#dest*time
model <- lme(grate ~ desttime, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(desttime="Tukey")), test = adjusted(type = "bonferroni"))

##UNTRANSFORMED WITH GROUPED ORIGIN_DEST
grate.aov<-aov(grate~ time*origin_dest + Error(colony), data=por)
summary(grate.aov)
#origin_dest
model <- lme(grate ~ origin_dest, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(origin_dest="Tukey")), test = adjusted(type = "bonferroni"))
#ORIGIN_dest*time
model <- lme(grate ~ grtime, random = ~ 1 | colony, por, na.action=na.exclude)
summary(glht(model, linfct=mcp(grtime="Tukey")), test = adjusted(type = "bonferroni"))

#####Goniastrea####
gon<-growth[growth$species=='gon',]
gon <- subset(gon, origin_dest != "LV_LV.Jan-16")#subsetting to remove since there is no data and messes up post hoc comparisons
#interactions for post hoc comparisons
gon$grtime <- factor(gon$grtime)
gon$desttime <- factor(gon$desttime)
summary(gon)

#looking for outliers
boxplot(grate ~ origin_dest + time, data=gon, las=2)
#shapiro
aggregate(grate ~ time + origin + dest, data=gon, FUN=function(x) shapiro.test(x)$p.value)
#bartletts HOV
for(i in c("16-Jan","16-Jul")){
	print(paste0("Timepoint",i))
	print(bartlett.test(grate ~ origin_dest, data=gon[gon$time==i,]))
}
plotNormalHistogram(gon$grate)
resid<-lm(grate~origin*dest*time, data=gon)
plot(residuals(resid)) #want random scatter; no apparent trendline
summary(resid)
qqPlot(resid,data=gon) 


##UNTRANSFORMED WITH ORIGIN DEST SEPARATE##
###repeated measures aov -> this is the same syntax for including colony as a source of error, would replace with time if running rm
grate.aov<-aov(grate ~ time*origin*dest + Error(colony), data=gon)
summary(grate.aov)
#dest
model <- lme(grate ~ dest, random = ~ 1 | colony, gon, na.action=na.exclude)
summary(glht(model, linfct=mcp(dest="Tukey")), test = adjusted(type = "bonferroni"))
#dest*time
model <- lme(grate ~ oritime, random = ~ 1 | colony, gon, na.action=na.exclude)
summary(glht(model, linfct=mcp(oritime="Tukey")), test = adjusted(type = "bonferroni"))

##UNTRANSFORMED WITH GROUPED ORIGIN_DEST##
grate.aov<-aov(grate ~ time*origin_dest + Error(colony), data=gon)
summary(grate.aov)
##To generate a priori contrasts in an unbalanced random mixed model
#origin_dest*time
miss1_mod<-lmer(grate ~ grtime + (1 | colony), gon, na.action=na.exclude)
anova(miss1_mod)
summary(glht(miss1_mod,linfct=mcp(grtime="Tukey")))
#origin_dest
miss1_mod<-lmer(grate ~ origin_dest + (1 | colony), gon, na.action=na.exclude)
anova(miss1_mod)
summary(glht(miss1_mod,linfct=mcp(origin_dest="Tukey")))

#####FIGURES####
#Rxn norm plot using ggplot
#POR#
sum=Summarize(grate~origin_dest+time, data=por, digits=3)
sum$se = sum$sd / sqrt(sum$n)
sum$se = signif(sum$se, digits=3)
sum$Origin<-factor(c("HV","MV","LV","MV","LV","HV","MV","LV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","LV","HV","HV","HV","MV","LV"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$dest=factor(sum$dest,levels=unique(sum$dest))
sum$time=factor(c("Jan-2016","Jan-2016","Jan-2016","Jan-2016","Jan-2016","Jul-2016","Jul-2016","Jul-2016","Jul-2016","Jul-2016"))
sum$time=factor(sum$time,levels=unique(sum$time))
sum

#poster and publication figure
library(cowplot)
pd=position_dodge(.75)
plot.growth<-ggplot(sum,aes(x=dest,y=mean,color=Origin)) + 
	geom_errorbar(aes(ymin=mean - se,ymax=mean + se), width=0.2,position=pd) + 	
  geom_point(aes(shape=Origin),size=4,position=pd) + 
	facet_wrap(~time) +
	theme_bw() + 
	theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'),
		panel.grid.minor = element_blank()) + 				 								
	scale_color_manual(values=c("red","gold","blue")) + 					
	coord_cartesian(xlim=c(0.5,3.25), ylim=c(0,0.4), expand=F) + 
	ylab("Weekly Growth (g/wk)") + 
	xlab("")
save_plot("20200417_por_growthsd.pdf", plot.growth,
          base_aspect_ratio = 1.3)

#GON#
#Interaction plot with transplant site on x axis.
sum=Summarize(grate~origin_dest+time, data=gon, digits=3)
sum$Origin<-factor(c("HV","MV","LV","MV","LV","HV","MV","LV","MV","LV"))
sum$dest<-factor(c("HV","HV","HV","MV","LV","HV","HV","HV","MV","LV"))
sum$Origin=factor(sum$Origin,levels=unique(sum$Origin))
sum$Time=factor(c("Jan-2016","Jan-2016","Jan-2016","Jan-2016","Jan-2016","Jul-2016","Jul-2016","Jul-2016","Jul-2016","Jul-2016"))
sum$dest<-factor(c("HV","HV","HV","MV","LV","HV","HV","HV","MV","LV"))
sum$dest=factor(sum$dest,levels=unique(sum$dest))
sum$time=factor(sum$time,levels=unique(sum$time))
sum$se = sum$sd / sqrt(sum$n)
sum$se = signif(sum$se, digits=3)
sum

#poster and publication figure
pd=position_dodge(.75)
plot.growth<-ggplot(sum,aes(x=dest,y=mean,color=Origin)) + 
	geom_errorbar(aes(ymin=mean - se,ymax=mean + se), width=0.2, size=0.7,position=pd) + 	geom_point(aes(shape=Origin),size=4,position=pd) + 
	facet_wrap(~Time) +
	theme_bw() + 
	theme(plot.margin=margin(0.5,0.5,0.5,0.5,'cm'),
		panel.grid.minor = element_blank()) + 				 								
		scale_colour_manual(values=c("red","gold","blue")) + 					
		coord_cartesian(xlim=c(0.5,3.3), ylim=c(0,0.15), expand=F) + 
	ylab("Weekly Growth (g/wk)") + 
	xlab("")
save_plot("20200417_gon_growthse_poster.pdf", plot.growth,
          base_aspect_ratio = 1.3)

####MISC####
#scatterplot of growth in the HV pool on y against growth in native pool on x axis
growth<-read.csv("fig.csv", header=T)
#set factor levels
growth$origin=factor(growth$origin,levels=unique(growth$origin))
growth$dest=factor(growth$dest,levels=unique(growth$dest))
growth$origin_dest=factor(growth$origin_dest,levels(growth$origin_dest)[c(1,3,2)])
growth$dest=factor(growth$dest,levels(growth$dest)[c(1,3,2)])
por<-growth[growth$species=='por',]
gon<-growth[growth$species=='gon',]
names(gon)
gon <- subset(gon, origin_dest != "LV_LV.Jan-16")#subsetting to remove since there is no data and messes up post hoc comparisons

por.scatter<-ggplot() +
  geom_point(por,mapping=aes(x=nativegrate,y=Hvgrate,color=origin,shape=time),size=4) +
  theme_bw() + 
  theme(plot.margin=margin(0.5,0,0,0.5,'cm'),
        panel.grid.minor = element_blank()) +
  geom_abline(slope=1,intercept=0) +
  scale_colour_manual(values=c("red","gold","blue")) + 	
  ylab("Growth in HV pool (g/wk)") + 
  xlab("Growth in native pool (g/wk)") 

gon.scatter<-ggplot() +
  geom_point(gon,mapping=aes(x=nativegrate,y=Hvgrate,color=origin,shape=time),size=4) +
  theme_bw() + 
  theme(plot.margin=margin(0.5,0,0,0.5,'cm'),
        panel.grid.minor = element_blank()) +
  geom_abline(slope=1,intercept=0) +
  scale_colour_manual(values=c("red","gold","blue")) + 	
  ylab(NULL) + 
  xlab("Growth in native pool (g/wk)") +
  theme(legend.position = "none")
#align yaxis scale
por.scatter = por.scatter + ylim(0, 0.45)
gon.scatter = gon.scatter + ylim(0, 0.45)
all.scatter<-plot_grid(por.scatter+theme(legend.position = 'none'), gon.scatter, align = 'vh', labels = c("A", "B"), hjust = -1)
#cowplot call legend from por.scatter
leg <- get_legend(por.scatter+theme(legend.box.background = element_rect(),
                                    legend.box.margin=margin(0,0,0,0,'cm'),
                                    legend.position= "bottom", 
                                    legend.box="horizontal" ,
                                    legend.text = element_text(size=9), 
                                    legend.title=element_text(size=10)) + labs(shape = "Time",color = "Origin"))
#put all together
#add the legend to a second row
all<-plot_grid(all.scatter, leg, nrow=2,rel_heights = c(6,1))
save_plot("20191023_growth-scatter.pdf", all)
