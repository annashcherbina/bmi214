library(ggplot2)
library(data.table)
library(reshape2)
data1=data.frame(read.table("sitesCa1.euc",sep='\t',header=F))
names(data1)=c("Step","BoundToCalcium","NotBoundToCalcium","MutatedSites")
data2=data.frame(read.table("sitesCa2.euc",sep='\t',header=F))
names(data2)=c("Step","BoundToCalcium","NotBoundToCalcium","MutatedSites")
m1=melt(data1,id="Step")
m2=melt(data2,id="Step")
ggplot(data=m1,aes(x=Step,y=value,group=variable,colour=variable))+
  geom_line(size=2)+
  xlab("Time Step")+
  ylab("Euclidean Distance")+
  ggtitle("Distance Atoms 17 - 96")+
  theme_bw(20)

ggplot(data=m2,aes(x=Step,y=value,group=variable,colour=variable))+
  geom_line(size=2)+
  xlab("Time Step")+
  ylab("Euclidean Distance")+
  ggtitle("Distance Atoms 298 - 385")+
  theme_bw(20)

#######QUESTIONS 22, 23 ######################
data3=data.frame(read.table("1FW4_med_out.features",header=T,sep='\t'))
data4=data.frame(read.table("1FW4_noCa_med_out.features",header=T,sep='\t'))
data5=data.frame(read.table("1FW4_hot_out.features",header=T,sep='\t'))
names(data3)=c("TimeStamp","Centers1_MED","Centers2_MED")
names(data4)=c("TimeStamp","Centers1_NoCa","Centers2_NoCa")
names(data5)=c("TimeStamp","Centers1_HOT","Centers2_HOT")
merged1=merge(data3,data4,by="TimeStamp")
m1=melt(merged1,id="TimeStamp")
merged2=merge(data3,data5,by="TimeStamp")
m2=melt(merged2,id="TimeStamp")
ggplot(data=m1,aes(x=TimeStamp,y=value,group=variable,colour=variable))+
  geom_line()+
  geom_point(size=2)+
  xlab("Time Step")+
  ylab("Score")+
  ggtitle("Comparison of Ca Binding Site Scores")+
  theme_bw(20)
ggplot(data=m2,aes(x=TimeStamp,y=value,group=variable,colour=variable))+
  geom_line()+
  geom_point(size=2)+
  xlab("Time Step")+
  ylab("Score")+
  ggtitle("Comparison of MED (350K) and HOT (4500K) Temperatures")+
  theme_bw(20)
