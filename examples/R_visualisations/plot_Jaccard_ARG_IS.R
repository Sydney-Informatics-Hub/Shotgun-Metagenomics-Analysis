
####################################################################################
## jaccard plot for Resistome, Insertion_sequences and Microbiome
####################################################################################
##### jaccard function
jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
}

##### plot

library(ggplot2)
library(reshape2)

dtp<-read.delim("jcc_arg_is_sp.txt")
d_tp<-melt(dtp, id.vars="timepoint")
names(d_tp) <- c("Timepoint", "variable","Jaccard_index")

d_tp<-d_tp[c(1:12,13:24,25:36),]
d_tp$type<- c(rep("Resistome",12),rep("Insertion_sequences",12),rep("Microbiome",12))

d_tp$groups<-c(rep("healthy",3),rep("mild",3),rep("moderate",3),rep("severe",3),rep("healthy",3),rep("mild",3),rep("moderate",3),rep("severe",3),rep("healthy",3),rep("mild",3),rep("moderate",3),rep("severe",3))


ggplot(d_tp, aes(Timepoint, Jaccard_index, group = variable, shape = type,color = groups,)) + geom_line(linetype = c(rep("dotted",12),rep("solid",12),rep("dotdash",12))) + theme_bw()+  geom_point(size=2.5)+scale_color_manual(values=c('dark green','dark orange','purple','red','black','purple','red','#E69F00','green','blue','black','purple','red','#E69F00','purple','blue','black','purple','red','#E69F00','green','blue','black','purple'))+theme(legend.text = element_text(size=8),legend.key.size = unit(0.5, 'cm'))




"#E64626", "#0148A4", "#FFB800", "#007E3B"


ggplot(d_tp, aes(Timepoint, Jaccard_index, group = variable, shape = type,color = groups,)) + geom_line(linetype = c(rep("dotted",12),rep("solid",12),rep("dotdash",12))) + theme_bw()+  geom_point(size=2.5)+scale_color_manual(values=c("#E64626", "#0148A4", "#FFB800", "#007E3B"))+theme(legend.text = element_text(size=8),legend.key.size = unit(0.5, 'cm'))