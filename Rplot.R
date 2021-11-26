# library
library(ggridges)
library(ggplot2)

# Diamonds dataset is provided by R natively
#head(diamonds)
#https://www.r-graph-gallery.com/294-basic-ridgeline-plot.html

filenames <- list.files(pattern="histogram_*", full.names=TRUE)

nd<-data.frame(name=NULL, val=NULL)
for (i in 1:length(filenames)){
  print(i)
  dd<-read.csv(filenames[i],sep=',', header = FALSE)
  dd=unlist(dd, use.names=FALSE)
  qd<-data.frame(name=rep(i,length(dd)), val=dd)
  nd<-rbind(nd,qd)
  
}

 
# basic example
ggplot(nd, aes(x = val, y = name, group=name, fill=name)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")

# basic example
ggplot(diamonds, aes(x = price, y = cut, fill = cut)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
