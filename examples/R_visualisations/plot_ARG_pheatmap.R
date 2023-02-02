###############################################################################################################
################## pheatmap of ARGs across samples of T1, T2, T3 (pres/abs)
###############################################################################################################

mat<-read.delim("ARG_pres_abs_sorted_530.txt") # 170 ARGs; pres.abs info

row.names(mat) <- mat$Sample_ID
mat<-mat[,4:173]

for (i in 1:length(mat)){
	mat[,i]<-as.numeric(mat[,i])
}

labs.row <- rownames(mat)

labs.row[1:530] <- ""
labs.row[69] <- "T1"
labs.row[228] <- "T2"
labs.row[425] <- "T3"

row<-read.delim("ARG_pres_abs_sorted_530.txt") # 170 ARGs; pres.abs info

df.row<-data.frame(row.names = row$Sample_ID, caries_degree =as.character(row$caries_states))

names(df.row)<-"caries_degree"
mycolors_row<-c( "#00CC99","#FFFF00","#CC0000", "#660000")
names(mycolors_row) <- c("Healthy","Mild","Moderate","Severe")

#### anno_row for ARG_class
class<-read.csv("ARG_class_renaming.csv")
m<-match(names(row)[4:173],class$ARG)

class_true<-read.csv("ARG_class.csv")
names(mat)<-class_true$ARG[m] ## to correct the special symbols in the ARG names

df.col<-data.frame(row.names = class_true$ARG[m], ARG_class= class_true$class[m])
names(df.col)<-"AMR_gene_class"

#library(randomcoloR)

n <- length(unique(class$class[m]))
palette <- distinctColorPalette(n)
mycolors_col<-palette 
# palette
# [1] "#DA634A" "#786FD8" "#CBDE77" "#D5AA96" "#CFE2DF" "#A1DEAD" "#D648DF" "#DABDDE"
# [9] "#70E4DB" "#D1A657" "#DB4CA5" "#84E950" "#8F6F94" "#E38BA5" "#7FA0E0" "#DC91E3"
#[17] "#7EBAD3" "#60E390" "#6F8972" "#E7E5B2" "#E2E33F" "#7B3EE1"

names(mycolors_col) <- unique(class$class[m])

ann_colors = list(
    caries_degree  = mycolors_row,
    AMR_gene_class = mycolors_col
)    


###############
#library(pheatmap)

## anno for both caries_degree and ARG_class
pheatmap(mat,gaps_row=c(139,319),annotation_row= df.row,annotation_col= df.col,cluster_rows = FALSE,fontsize_row = 12, fontsize_col = 3, fontsize=7, annotation_colors= ann_colors,labels_row = labs.row,clustering_method = "ward.D2") # "ward.D2"



