library(ggridges)
library(ggplot2)
library(dplyr)
library(reshape2)
library(moments)
library(ggpubr)

WD <- getwd()

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

data_process<-function(location, name){
  nd<-data.frame(name=NULL, val=NULL)
  nd2<-NULL
  for (j in 1:length(location)){
    setwd(location[j])
    filenames <- list.files(pattern="^histogram_pd_*", full.names=TRUE)
    filenames2 <- list.files(pattern=paste0('init_core_*'), full.names=TRUE)
    
    for (i in 1:length(filenames)){
      print(i)
      cond=name[j]
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
      
      core=gsub('.*_([0-9]+).*','\\1', filenames[i])
      qd<-data.frame(name=rep(core,length(dD)), D=dD, mD=dmD, G3P=dg3p, mR=dmR, R=dR, pD=dpD, pR=dpR, cycle=as.character(mean(dd2$V1)), cond=cond)
      nd<-rbind(nd,qd)
      
      mean_d<-data.frame(name=core, D=mean(dD), R=mean(dR), mD=mean(dmD), mR=mean(dmR), 
                         cv_D=sd(dD)/mean(dD), cv_R=sd(dR)/mean(dR), cv_mD=sd(dmD)/mean(dmD), cv_mR=sd(dmR)/mean(dmR), 
                         cycle=as.character(mean(dd2$V1)),cond=cond)
      
      nd2<-rbind(nd2,mean_d)
      
    }
    setwd(WD)
  }
  
  my_list <- list("nd" = nd, "nd2" = nd2)
  return(my_list)
  
  
}


loc1<-"noGrowth_Memory_noSwitching/low_ino"
loc2<-"noGrowth_Memory_noSwitching/high_ino"

loc3<-"noGrowth_Memory_Switching/low_ino"
loc4<-"noGrowth_Memory_Switching/high_ino"

loc5<-"noGrowth_noMemory_Switching/low_ino"
loc6<-"noGrowth_noMemory_Switching/high_ino"

my_comparisons <- list(c('Low(M+NS) inoculum','High(M+NS) inoculum'),c('Low(M+SW) inoculum','High(M+SW) inoculum'),
                       c('Low(NM+SW) inoculum','High(NM+SW) inoculum'),
                       c('Low(M+NS) inoculum','Low(M+SW) inoculum'),
                       c('Low(M+NS) inoculum','Low(NM+SW) inoculum'),
                       c('High(M+SW) inoculum','High(NM+SW) inoculum'),
                       c('High(M+NS) inoculum','High(M+SW) inoculum'),
                       c('High(M+NS) inoculum','High(NM+SW) inoculum'),
                       c('High(M+SW) inoculum','High(NM+SW) inoculum'))


loc<-c(loc1, loc2, loc3, loc4, loc5, loc6)
name<-c('Low(M+NS) inoculum','High(M+NS) inoculum','Low(M+SW) inoculum','High(M+SW) inoculum','Low(NM+SW) inoculum','High(NM+SW) inoculum')#,'Low inoculum and growth modulation','High inoculum and growth modulation')
d1<-data_process(loc,name)

nd<-d1$nd
nd2<-d1$nd2

mc_10_m<-nd %>% group_by(name, cond) %>%   summarise(mean = mean(mD), sd = sd(mD), cv=sd(mD)/mean(mD) )
den_d<-melt(nd[,c(1,3,10)])


#basic example
pdf('sample_plot.pdf', height=10, width = 18)
# png('glpd_dist.png')
p1<-ggplot(den_d[den_d$cond%in%c('Low(M+NS) inoculum','High(M+NS) inoculum'),], aes(x = value, y = name,group=name, fill=name)) +
  geom_density_ridges(bandwidth = 0.5,quantile_lines = TRUE, quantiles = 4,scale=3,alpha = .7,vline_color=alpha('black',1)) +
  theme_ridges()+facet_wrap(~cond,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold", size=16,angle=45, vjust = 0.5),
    legend.position = "none",
    axis.text.y = element_text(face="bold", size=16),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 16, color = "black", face = "bold.italic"))+xlab('GlpD')+xlim(c(-1,25))
#  xlim(c(-1,30))+ylab('Colony')

p2<-ggplot(den_d[den_d$cond%in%c('Low(M+SW) inoculum','High(M+SW) inoculum'),], aes(x = value, y = name,group=name, fill=name)) +
  geom_density_ridges(bandwidth = 0.5,quantile_lines = TRUE, quantiles = 4,scale=3,alpha = .7,vline_color=alpha('black',1)) +
  theme_ridges()+facet_wrap(~cond,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold", size=16,angle=45, vjust = 0.5),
    legend.position = "none",
    axis.text.y = element_text(face="bold", size=16),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 16, color = "black", face = "bold.italic"))+xlab('GlpD')+xlim(c(-1,25))

p3<-ggplot(den_d[den_d$cond%in%c('Low(NM+SW) inoculum','High(NM+SW) inoculum'),], aes(x = value, y = name,group=name, fill=name)) +
  geom_density_ridges(bandwidth = 0.5,quantile_lines = TRUE, quantiles = 4,scale=3,alpha = .7,vline_color=alpha('black',1)) +
  theme_ridges()+facet_wrap(~cond,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold", size=16,angle=45, vjust = 0.5),
    legend.position = "none",
    axis.text.y = element_text(face="bold", size=16),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 16, color = "black", face = "bold.italic"))+xlab('GlpD')+xlim(c(-1,25))


#mylayout <- matrix(seq(1,2,3), cols = 3)

multiplot(p1,p2,p3, cols =3 )

dev.off()

# 
# colnames(nd2)<-c("name",  "mean D",     "mean R",     "mean mD",    "mean mR",    "cv D",  "cv R",  "cv mD", "cv mR", "cycle", "cond") 
# p_p<-melt(nd2)
# new_dat<-p_p[which(p_p$variable %in% c('mean D', 'mean mD','cv D','cv mD')),]
# 
# pdf('stat_tt2.pdf')
# ggplot(new_dat, aes(x=cond,y=value)) + 
#   geom_boxplot()+#coord_trans(y="log10")+
#   geom_jitter(shape=16,col='red',position=position_jitter(0.1))+facet_wrap(~variable,scales="free_y")+
#   theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
#     panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#     legend.position = c(0.8, 0.2),legend.key.size = unit(1, 'cm'),legend.title = element_blank(),
#     legend.text = element_text(colour="black", size=15), axis.text.x = element_text(face="bold", size=12,angle=45, vjust = 0.5),
#     axis.text.y = element_text(face="bold", size=12),
#     axis.title.x = element_text(face="bold", size=16),
#     axis.title.y = element_text(face="bold", size=16),
#     strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
#     strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+xlab('Condition')+ylab('')
# dev.off()

kur<-nd %>% group_by(name, cond) %>%   summarise('Mean glpD'=mean(mD), CV=sd(mD)/mean(mD), Skewness = skewness(mD), Kurtosis = kurtosis(mD) )
kur<-melt(kur, id=c('name','cond'))
pdf('stat_sample.pdf', height = 7)
ggplot(kur, aes(x=cond,y=value)) + 
  geom_boxplot()+#coord_trans(y="log10")+
  geom_jitter(shape=16,col='red',position=position_jitter(0.1))+
  stat_compare_means(comparisons = my_comparisons)+
  facet_wrap(~variable,scales="free_y")+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = c(0.8, 0.2),legend.key.size = unit(1, 'cm'),legend.title = element_blank(),
    legend.text = element_text(colour="black", size=15), axis.text.x = element_text(face="bold", size=10,angle=90, vjust = 0.5),
    axis.text.y = element_text(face="bold", size=12),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+xlab('')+ylab('')+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 14))
dev.off()
