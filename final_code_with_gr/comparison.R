#latest R file
# library
library(ggridges)
library(ggplot2)
library(dplyr)
library(reshape2)
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# Diamonds dataset is provided by R natively
#head(diamonds)
#https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html
#
#filenames <- list.files(path="good_result/king3p/low_copy/",pattern="^histogram_pd_*", full.names=TRUE)
#filenames2 <- list.files(path="good_result/king3p/low_copy/",pattern=paste0('init_core_*'), full.names=TRUE)

loc1<-'/home/abp19/Projects/transcriptional_noise/good_result/Sus_model/inves/final/testing2/exp0.05/10/'
loc2<-'/home/abp19/Projects/transcriptional_noise/good_result/Sus_model/inves/final/testing2/exp0.05/100/'
setwd(loc1)
filenames <- list.files(pattern="^histogram_pd_*", full.names=TRUE)
#filenames <- list.files(path="simple_model/low_copy/indep_5/", pattern="^histogram_pd_*", full.names=TRUE)
#filenames<-filenames[c(1:10)]
#filenames <- list.files(pattern=paste0('_histogram_100.0_*'), full.names=TRUE)
filenames2 <- list.files(pattern=paste0('init_core_*'), full.names=TRUE)

#filenames <- paste0('./histogram_pd_',seq(1:20),'.csv')
#filenames2 <- paste0('./init_core_',seq(1:20),'_histogram_0.1.csv')

nd<-data.frame(name=NULL, val=NULL)
nd2<-NULL
for (i in 1:length(filenames)){
  print(i)
  cond='low'
  dd<-read.csv(filenames[i],sep=',', header = FALSE)
  # dd<-dd[sample(nrow(dd),1000),]
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
  qd<-data.frame(name=rep(core,length(dD)), D=dD, mD=dmD, G3P=dg3p, mR=dmR, R=dR, pD=dpD, pR=dpR, cycle=as.character(mean(dd2$V1)), cond=cond)
  nd<-rbind(nd,qd)
  # mean_d<-data.frame(name=core, D=median(dD), mD=median(dmD), G3P=median(dg3p), mR=median(dmR), R=median(dR), pD=median(dpD), pR=median(dpR), 
  #            sd_D=sd(dD), sd_mD=sd(dmD), sd_G3P=sd(dg3p), sd_mR=sd(dmR), sd_R=sd(dR), sd_pD=sd(dpD), sd_pR=sd(dpR),
  #            cv_D=sd(dD)/mean(dD), cv_mD=sd(dmD)/mean(dmD), cv_G3P=sd(dg3p)/mean(dg3p), cv_mR=sd(dmR)/mean(dmR), cv_R=sd(dR)/mean(dR), 
  #            cv_pD=sd(dpD)/mean(dpD), cv_pR=sd(dpR)/mean(dpR),cycle=as.character(mean(dd2$V1)))
  
  mean_d<-data.frame(name=core, D=mean(dD), R=mean(dR), mD=mean(dmD), mR=mean(dmR), 
                     cv_D=sd(dD)/mean(dD), cv_R=sd(dR)/mean(dR), cv_mD=sd(dmD)/mean(dmD), cv_mR=sd(dmR)/mean(dmR), 
                     cycle=as.character(mean(dd2$V1)),cond=cond)
  
  nd2<-rbind(nd2,mean_d)
  
}

setwd(loc2)

filenames <- list.files(pattern="^histogram_pd_*", full.names=TRUE)
filenames2 <- list.files(pattern=paste0('init_core_*'), full.names=TRUE)
for (i in 1:length(filenames)){
  print(i)
  cond='high'
  dd<-read.csv(filenames[i],sep=',', header = FALSE)
  # dd<-dd[sample(nrow(dd),1000),]
  
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
  qd<-data.frame(name=rep(core,length(dD)), D=dD, mD=dmD, G3P=dg3p, mR=dmR, R=dR, pD=dpD, pR=dpR, cycle=as.character(mean(dd2$V1)), cond=cond)
  nd<-rbind(nd,qd)
  # mean_d<-data.frame(name=core, D=median(dD), mD=median(dmD), G3P=median(dg3p), mR=median(dmR), R=median(dR), pD=median(dpD), pR=median(dpR), 
  #            sd_D=sd(dD), sd_mD=sd(dmD), sd_G3P=sd(dg3p), sd_mR=sd(dmR), sd_R=sd(dR), sd_pD=sd(dpD), sd_pR=sd(dpR),
  #            cv_D=sd(dD)/mean(dD), cv_mD=sd(dmD)/mean(dmD), cv_G3P=sd(dg3p)/mean(dg3p), cv_mR=sd(dmR)/mean(dmR), cv_R=sd(dR)/mean(dR), 
  #            cv_pD=sd(dpD)/mean(dpD), cv_pR=sd(dpR)/mean(dpR),cycle=as.character(mean(dd2$V1)))
  
  mean_d<-data.frame(name=core, D=mean(dD), R=mean(dR), mD=mean(dmD), mR=mean(dmR), 
                     cv_D=sd(dD)/mean(dD), cv_R=sd(dR)/mean(dR), cv_mD=sd(dmD)/mean(dmD), cv_mR=sd(dmR)/mean(dmR), 
                     cycle=as.character(mean(dd2$V1)),cond=cond)
  
  nd2<-rbind(nd2,mean_d)
  
}
# ggplot(nd, aes(x = D, y = cycle, group=cycle, fill=cycle)) +
#   geom_density_ridges() +
#   theme_ridges()+
#   theme(legend.position = "none")

#pdf('10_gr_30_90.pdf',height=15, width = 8)
#ggplot(nd, aes(x =mD, y = name, group=name, fill=name)) +
#  geom_density_ridges(bandwidth = 0.5,alpha=.7) +
#  theme_ridges()+ 
#  theme(legend.position = "none")

#ggplot(nd)+ geom_density_ridges(aes(x = mD, y = name, group=name, fill=name,height = after_stat(count)), stat="density",
#                                scale = 1.5,
#                                alpha = 0.5) +
#  theme_ridges()+ 
#  theme(legend.position = "none")
#dev.off()
#dev.new()
#ggplot(nd, aes(x = mD, y = name, group=name, fill=name,height = ..density..)) +
#  geom_density_ridges(alpha=.7,stat = "density") +
#  theme_ridges()+ 
#  theme(legend.position = "none")
#ggplot(nd, aes(x = mD, y = name, group=name, fill=name)) +
#  geom_density_ridges(stat='binline',alpha=.7) +
#  theme_ridges()+ 
#  theme(legend.position = "none")

mc_10_m<-nd %>% group_by(name, cond) %>%   summarise(mean = mean(mD), sd = sd(mD), cv=sd(mD)/mean(mD) )
den_d<-melt(nd[,c(1,3,10)])

setwd('/home/abp19/Projects/transcriptional_noise/good_result/Sus_model/inves/final/')

#basic example
pdf('comparison_testing.pdf', height=15, width = 12)
# png('glpd_dist.png')
ggplot(den_d, aes(x = value, y = name,group=name, fill=name)) +
  geom_density_ridges(bandwidth = 0.5,quantile_lines = TRUE, quantiles = 4,scale=3,alpha = .7,vline_color=alpha('black',1)) +
  theme_ridges()+facet_wrap(~cond,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
     axis.text.x = element_text(face="bold", size=12,angle=45, vjust = 0.5),
    legend.position = "none",
    axis.text.y = element_text(face="bold", size=12),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+xlab('GlpD')#+
#  xlim(c(-1,30))+ylab('Colony')

#ggplot(den_d, aes(x = value, y = name, fill=cond)) +
#  geom_density_ridges(scale=3,alpha = .7,vline_color=alpha('black',0.7)) +
#  theme_ridges()+facet_wrap(~cond,scales="free_y")+
#  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
#    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
#    axis.text.x = element_text(face="bold", size=12,angle=45, vjust = 0.5),
#    legend.position = 'top',legend.text = element_text(colour="black", size=15),legend.title=element_blank(),
#    plot.title = element_text(color="red", size=14, face="bold.italic"),
#    axis.text.y = element_text(face="bold", size=12),
#    axis.title.x = element_text(face="bold", size=16),
#    axis.title.y = element_text(face="bold", size=16),
#    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
#    strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+xlab('GlpD')+xlim(c(-5,75))+ylab('Colony')+ggtitle('')
dev.off()


colnames(nd2)<-c("name",  "mean D",     "mean R",     "mean mD",    "mean mR",    "cv D",  "cv R",  "cv mD", "cv mR", "cycle", "cond") 
p_p<-melt(nd2)
new_dat<-p_p[which(p_p$variable %in% c('mean D', 'mean mD','cv D','cv mD')),]

pdf('stat_testing.pdf')
ggplot(new_dat, aes(x=cond,y=value)) + 
  geom_boxplot()+#coord_trans(y="log10")+
  geom_jitter(shape=16,col='red',position=position_jitter(0.1))+facet_wrap(~variable,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = c(0.8, 0.2),legend.key.size = unit(1, 'cm'),legend.title = element_blank(),
    legend.text = element_text(colour="black", size=15), axis.text.x = element_text(face="bold", size=12,angle=45, vjust = 0.5),
    axis.text.y = element_text(face="bold", size=12),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+xlab('Condition')+ylab('')
dev.off()






stop()

ggplot(den_d)+geom_density_ridges(aes(x = value, y = name,group=name, fill=name,height = after_stat(count)),stat="density",
                                  scale = 1.5,
                                  alpha = 0.5,vline_color=alpha('black',1))+
  theme_ridges()+facet_wrap(~cond,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold", size=12,angle=45, vjust = 0.5),
    legend.position = "none",
    axis.text.y = element_text(face="bold", size=12),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+xlab('GlpD')#+
xlim(c(-1,30))+ylab('Colony')


p2<-ggplot(den_d, aes(x = value, y = name, fill=cond)) +
  geom_density_ridges(scale=3,alpha = .7,vline_color=alpha('black',0.7)) +
  theme_ridges()+#facet_wrap(~cond,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold", size=12,angle=45, vjust = 0.5),
    legend.position = 'top',legend.text = element_text(colour="black", size=15),legend.title=element_blank(),
    plot.title = element_text(color="red", size=14, face="bold.italic"),
    axis.text.y = element_text(face="bold", size=12),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+xlab('GlpD')+
  xlim(c(-5,75))+ylab('Colony')+ggtitle('')


#den_d_h<-den_d

mylayout <- matrix(seq(1,2), ncol = 2, byrow = TRUE)

pdf('add_g3p_low_high.pdf', width = 12, height = 16)
multiplot(p1,p2, layout = mylayout )
dev.off()






# filenames <- list.files(path='low_copy/',pattern="histogram_pd_*", full.names=TRUE)
# 
# filenames2 <- list.files(pattern=paste0('_histogram_300.0_*'), full.names=TRUE)
# nd<-data.frame(name=NULL, val=NULL)
# cv<-NULL
# for (i in 1:length(filenames2)){
#   print(i)
#   # filenames <- list.files(pattern=paste0("tnow_core_",i,'_histogram_50.0_*'), full.names=TRUE)
#   filenames<-filenames2[i]
#   nd2<-NULL
#   for (j in 1:length(filenames)){
#     dd<-read.csv(filenames[j],sep=',', header = FALSE)
#     nd2<-rbind(nd2,dd)
#   }
# 
#   qd<-data.frame(name=rep(i,length(nd2$V1)), val=nd2$V1)
#   print(length(qd$val))
#   nd<-rbind(nd,qd)
#   cv=rbind(cv,c(qd[1,1],mean(qd$val),sd(qd$val),sd(qd$val)/mean(qd$val)))
# }
# cv<-data.frame(cv)
# colnames(cv)<-c('Name', 'mean', 'sd', 'cv')
# 
# # basic example
# ggplot(nd, aes(x = val, y = name, group=name, fill=name)) +
#   geom_density_ridges() +
#   theme_ridges()+
#   theme(legend.position = "none")
# # p1<-p1+xlim(c(0,60))
# p2<-p2+xlim(c(0,60))
# 
# library(ggpubr)
# ggarrange(p1,p2, ncol = 2, nrow = 1)
# 
# 
# 
# lc_dat<-cv
# lc_dat$name<-rep('low',length(lc_dat$mean))
# 
# hc_dat<-cv
# hc_dat$name<-rep('high',length(hc_dat$mean))
# 
# fd<-rbind(lc_dat,hc_dat)
# 
# ggplot(fd,aes(x=cv, colour=name)) + geom_density()


