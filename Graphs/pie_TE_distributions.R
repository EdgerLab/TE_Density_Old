# By Scott Teresi
# Sample Pie Graph
#=============================================
library(ggplot2)
library(RColorBrewer)
library(plotrix)
setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/")
#Type_TE_Data = read.csv2('Type_Table.csv', header=FALSE)
#Family_TE_Data = read.csv2('Family_Table.csv', header=FALSE)

Type_TE_Data <- data.frame(
  group = c("DNA", "LTR", "Unknown", "LINE"),
  value = c(82619, 314043, 45600, 15646)
)
head(Type_TE_Data)

#colnames(Type_TE_Data) = c("Type","Number")
#colnames(Family_TE_Data) = c("Family","Number")

#bp<- ggplot(Type_TE_Data, aes(x="", y=value, fill=group))+
#geom_bar(width = 1, stat = "identity") # this step is required to get the colors and such

#pie <- bp + coord_polar("y", start=0)
#pie










slices = c(82619, 314043, 45600, 15646)
lbls = c("DNA", "LTR", "Unknown", "LINE")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="TE Percentages Genome-Wide by Type")


png(file="Type_Pie.png",width=700,height=600)
#plot(x=rnorm(10),y=rnorm(10),main="example")
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="TE Percentages Genome-Wide by Type")
par(mar=c(4,4,4,6)) # bottom, left top right
dev.off()



lbls = c("MULE","Gypsy (LTR)","Unknown LTR Family","CMC-EnSpm","Copia (LTR)","Unknown","LINE","hAT","PIF-Harbinger")
slices = c(17424,98511,174430,17233,41100,65293,15646,18657,9609)
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels
lbls <- paste(lbls,"%",sep="") # ad % to labels
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="TE Percentages Genome-Wide by Family")


png(file="Family_Pie.png",width=700,height=600)
#plot(x=rnorm(10),y=rnorm(10),main="example")
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="TE Percentages Genome-Wide by Family")
par(mar=c(4,4,4,6))
dev.off()

