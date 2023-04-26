# library
library(ggridges)
library(ggplot2)
library(dplyr)
# Diamonds dataset is provided by R natively
#head(diamonds)
#https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html

filenames <- list.files(path="simple_model/low_copy/gillespie_nparam_nmodel_3/",pattern="^histogram_pd_*", full.names=TRUE)

filenames <- list.files(pattern="^histogram_pd_*", full.names=TRUE)
#filenames <- list.files(path="simple_model/low_copy/indep_5/", pattern="^histogram_pd_*", full.names=TRUE)
#filenames<-filenames[c(1:10)]
#filenames <- list.files(pattern=paste0('_histogram_100.0_*'), full.names=TRUE)
filenames2 <- list.files(pattern=paste0('init_core_*'), full.names=TRUE)

nd<-data.frame(name=NULL, val=NULL)
nd2<-NULL
for (i in 1:length(filenames)){
  print(i)
  dd<-read.csv(filenames[i],sep=',', header = FALSE)
  print(length(dd$V1))
  dD=unlist(dd$V1, use.names=FALSE)
  dmD=unlist(dd$V2, use.names=FALSE)
  dg3p=unlist(dd$V3, use.names=FALSE)
  dmR=unlist(dd$V4, use.names=FALSE)
  dR=unlist(dd$V5, use.names=FALSE)
  dpD=unlist(dd$V6, use.names=FALSE)
  dpR=unlist(dd$V7, use.names=FALSE)
  
  dd2<-read.csv(filenames2[i],sep=',', header = FALSE)
  #print(c(mean(dd),sd(dd)))
  # nd2<-rbind(nd2,c(mean(dd),sd(dd),sd(dd)/mean(dd)))
  # core=i
  core=gsub('.*_([0-9]+).*','\\1', filenames[i])
  qd<-data.frame(name=rep(core,length(dD)), D=dD, mD=dmD, G3P=dg3p, mR=dmR, R=dR, pD=dpD, pR=dpR, cycle=as.character(mean(dd2$V1)))
  nd<-rbind(nd,qd)
  
}

ggplot(nd, aes(x = D, y = cycle, group=cycle, fill=cycle)) +
  geom_density_ridges() +
  theme_ridges()+
  theme(legend.position = "none")

ggplot(nd, aes(x = D, y = name, group=name, fill=name)) +
  geom_density_ridges() +
  theme_ridges()+ 
  theme(legend.position = "none")

mc_10_m<-nd %>% group_by(name) %>%   summarise(mean = mean(D), sd = sd(D), cv=sd(D)/mean(D) )

#basic example
pdf('t_100_m.pdf', height=15, width = 5)
ggplot(nd, aes(x = D, y = name, group=name, fill=name)) +
  geom_density_ridges() +
  theme_ridges()+ 
  theme(legend.position = "none")+xlim(c(-10,80))
dev.off()



lc_dat<-data.frame(nd2)

colnames(lc_dat)<-c('mean','sd','cv')
lc_dat$name<-rep('low',length(lc_dat$mean))

hc_dat<-data.frame(nd2)
colnames(hc_dat)<-c('mean','sd','cv')
hc_dat$name<-rep('high',length(hc_dat$mean))

fd<-rbind(lc_dat,hc_dat)

#

ggplot(nd, aes(x = val,y='tt')) +
  geom_density_ridges() +
  theme_ridges()+ 
  theme(legend.position = "none")

ggsave("high_copy_mrna_distribution.pdf", width = 4, height = 4)

ggplot(fd,aes(x=mean, colour=name)) + geom_density() 

ggplot(fd, aes(x=name, y=cv)) + 
  geom_boxplot()+geom_jitter(shape=16,col='red',position=position_jitter(0.1))

+ facet_wrap(~name)+xlim(0,150)

ld<-data.frame(dat=c(lc_dat$sd,hc_dat$sd),col=rep(c('low','high'),c(30,30)))

# basic example
ggplot(ld, aes(x = dat, color=col)) +
  geom_density()




filenames <- list.files(path='low_copy/',pattern="histogram_pd_*", full.names=TRUE)

filenames2 <- list.files(pattern=paste0('_histogram_150.0_*'), full.names=TRUE)
nd<-data.frame(name=NULL, val=NULL)
cv<-NULL
for (i in 1:length(filenames2)){
  print(i)
  # filenames <- list.files(pattern=paste0("tnow_core_",i,'_histogram_50.0_*'), full.names=TRUE)
  filenames<-filenames2[i]
  nd2<-NULL
  for (j in 1:length(filenames)){
    dd<-read.csv(filenames[j],sep=',', header = FALSE)
    nd2<-rbind(nd2,dd)
  }
  
  qd<-data.frame(name=rep(i,length(nd2$V1)), val=nd2$V1)
  print(length(qd$val))
  nd<-rbind(nd,qd)
  cv=rbind(cv,c(qd[1,1],mean(qd$val),sd(qd$val),sd(qd$val)/mean(qd$val)))
}
cv<-data.frame(cv)
colnames(cv)<-c('Name', 'mean', 'sd', 'cv')

# basic example
ggplot(nd, aes(x = val, y = name, group=name, fill=name)) +
  geom_density_ridges() +
  theme_ridges()+ 
  theme(legend.position = "none")
p1<-p1+xlim(c(0,60))
p2<-p2+xlim(c(0,60))

library(ggpubr)
ggarrange(p1,p2, ncol = 2, nrow = 1)



lc_dat<-cv
lc_dat$name<-rep('low',length(lc_dat$mean))

hc_dat<-cv
hc_dat$name<-rep('high',length(hc_dat$mean))

fd<-rbind(lc_dat,hc_dat)

ggplot(fd,aes(x=cv, colour=name)) + geom_density()


