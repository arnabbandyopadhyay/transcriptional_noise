# library
library(ggridges)
library(ggplot2)

# Diamonds dataset is provided by R natively
#head(diamonds)
#https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html

filenames <- list.files(path='low_copy/',pattern="histogram_pd_*", full.names=TRUE)
#filenames <- list.files(pattern="histogram_pd_*", full.names=TRUE)

nd<-data.frame(name=NULL, val=NULL)
nd2<-NULL
for (i in 1:length(filenames)){
  print(i)
  dd<-read.csv(filenames[i],sep=',', header = FALSE)
  print(length(dd$V1))
  dd=unlist(dd$V1, use.names=FALSE)
  #print(c(mean(dd),sd(dd)))
  nd2<-rbind(nd2,c(mean(dd),sd(dd),sd(dd)/mean(dd)))
  qd<-data.frame(name=rep(i,length(dd)), val=dd)
  nd<-rbind(nd,qd)
  
}


lc_dat<-data.frame(nd2)
colnames(lc_dat)<-c('mean','sd','cv')
lc_dat$name<-rep('low',length(lc_dat$mean))

hc_dat<-data.frame(nd2)
colnames(hc_dat)<-c('mean','sd','cv')
hc_dat$name<-rep('high',length(hc_dat$mean))

fd<-rbind(lc_dat,hc_dat)

# basic example
ggplot(nd, aes(x = val, y = name, group=name, fill=name)) +
  geom_density_ridges() +
  theme_ridges()+ 
  theme(legend.position = "none")

ggsave("high_copy_mrna_distribution.pdf", width = 4, height = 4)

ggplot(fd,aes(x=cv, colour=name)) + geom_density() 

+ facet_wrap(~name)+xlim(0,150)

ld<-data.frame(dat=c(lc_dat$sd,hc_dat$sd),col=rep(c('low','high'),c(30,30)))

# basic example
ggplot(ld, aes(x = dat, color=col)) +
  geom_density()
