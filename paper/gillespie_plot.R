library(ggplot2)
library(reshape2)

data<-read.csv('test.csv', header=T)[,c(1,2,5)]
data<-data[data$param %in% c(0.1,0.2,0.3,0.4,0.5, 0.7),]
data<-melt(data,id='param')
data$param<-as.factor(data$param)

pdf('gillespie.pdf',h=10,w=16)
ggplot(data, aes(x = param, y = value,color=variable,fill=variable)) +
  geom_violin(scale='width',alpha=0.5,adjust = 1.2)+
  annotate("text", label = c("Glucose",'LB','Glycerol'),x = c(0.8,3.5,6), y = c(16.5), size = 6,colour = "red",fontface = "bold.italic")+
  annotate("rect", xmin=c(0,1.5,5.5), xmax=c(1.5,5.5,6.5), ymin=c(0,0,0) , ymax=c(16), alpha=0.3, color="black", fill=c('grey','white','grey'))+
  geom_violin(scale='width',alpha=0.5,adjust = 1.2)+
  geom_boxplot( fill='white',position = position_dodge(0.9),width=0.1,show.legend=FALSE)+
  geom_boxplot( col='black',position = position_dodge(0.9),width=0.1,alpha=.001,show.legend=FALSE)+
  theme_bw()+theme(strip.background=element_rect(fill='white'))+theme(
    panel.background = element_rect(fill = "white",colour = "black",linewidth = 1),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    axis.text.x = element_text(face="bold", size=12,angle=45, vjust = 0.5),
    legend.position = c(0.1,0.8),legend.background = element_blank(),legend.title = element_blank(),
    legend.text = element_text( size=30, face="bold"),
    axis.text.y = element_text(face="bold", size=12),
    axis.title.x = element_text(face="bold", size=16),
    axis.title.y = element_text(face="bold", size=16),
    strip.text.x = element_text(size = 14, color = "black", face = "bold.italic"),
    strip.text.y = element_text(size = 14, color = "black", face = "bold.italic"))+ylab('glpD')+ylim(c(0,16.5))+ylab('mRNA')
#  xlim(c(-1,30))+ylab('Colony')
dev.off()