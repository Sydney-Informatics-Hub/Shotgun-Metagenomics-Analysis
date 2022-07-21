
###############################################################################################################
###### horizontal stacked bar for ACE model results
###############################################################################################################

library(ggplot2)
## prepare dataframe
table <- read.csv("results_ACE.csv")

data1 <- table[1:21,]
data2 <- table[24:44,]
data3 <- table[47:67,]

## generate stacked bar chart  
library(ggpubr)

p1 <- ggplot(data1, aes(fill=factor(ACE_category, levels=c("E","C","A")), y=as.numeric(value1), x= ARG)) + geom_bar(position="stack", stat="identity",width=0.5)+theme(axis.text=element_text(size=6,angle=90),legend.text = element_text(size=7),legend.key.size=unit(0.37, "cm"),legend.key.height=unit(0.4, "cm"),legend.position = c(1.05, 0.87))+scale_fill_manual("ACE_category", values =  c( "orange","dark green","purple"))+coord_flip()+theme_bw()+labs(y="")+labs(x="")

p2 <- ggplot(data2, aes(fill=factor(ACE_category, levels=c("E","C","A")), y=as.numeric(value2), x= ARG)) + geom_bar(position="stack", stat="identity",width=0.5)+theme(axis.text=element_text(size=6,angle=90),legend.text = element_text(size=7),legend.key.size=unit(0.37, "cm"),legend.key.height=unit(0.4, "cm"),legend.position = c(1.05, 0.87))+scale_fill_manual("ACE_category", values = c( "orange","dark green","purple"))+coord_flip()+theme_bw()+labs(y="")

p3 <- ggplot(data3, aes(fill=factor(ACE_category, levels=c("E","C","A")), y=as.numeric(value3), x= ARG)) + geom_bar(position="stack", stat="identity",width=0.5)+theme(axis.text.y=element_text(size=6,angle=30),legend.text = element_text(size=7),legend.key.size=unit(0.37, "cm"),legend.key.height=unit(0.4, "cm"),legend.position = c(1.05, 0.87))+scale_fill_manual("ACE_category", values = c( "orange","dark green","purple"))+coord_flip()+theme_bw()+labs(y="values")+labs(x="")

ggarrange(p1,p2,p3,ncol=1,nrow=3,common.legend=T, legend="right", align = "hv", widths = c(1, 1, 1))


