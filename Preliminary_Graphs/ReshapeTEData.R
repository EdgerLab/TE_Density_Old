#install.packages("reshape2")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("cowplot")
library("stringi")
library("tidyr")
library("dplyr")
library("reshape2")
library("cowplot")
#load Data
#`500_cleaned_data` <- read.csv("C:/Users/Kevin/Downloads/500_cleaned_data (2).csv", quote="'")
`500_cleaned_data` <- read.csv("500_cleaned_data.csv", quote="'")
#Reeshape data to make column for variable TE type/family and corresponding value for TE density
Reshaped<-melt(data = `500_cleaned_data`[-1],id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length"))
#Calculate averages by TE type/family
TE.Density.Means<-Reshaped %>% group_by(variable) %>% summarise(avg=mean(value)) %>% arrange(avg)

#Split TE type/family and location (e.g. left, right, intragenic) in hacky way (reverse string, split by "_", reverse strings again)
TE.Density.Means$variable <-stri_reverse(TE.Density.Means$variable)
TE.Density.Means<-TE.Density.Means%>%separate(variable,into = c("Location","TE_type"),sep="_",extra="merge")
TE.Density.Means$Location<-stri_reverse(TE.Density.Means$Location)
TE.Density.Means$TE_type<-stri_reverse(TE.Density.Means$TE_type)
#add Column for window size
TE.Density.Means$Distance<-"500bp"
#Make colourblind friendly palette :)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
myColors <- colorRampPalette(cbPalette)(12)
names(myColors) <- unique(TE.Density.Means$TE_type)
colScale <- scale_fill_manual(name="", values=myColors)
#Reorder Location factors so that plot shows left, intra, right
TE.Density.Means$Location <- factor(TE.Density.Means$Location, levels=c("left", "intra", "right"))
# Subset to look only at TE types
TE.Type.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="DNA"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="LTR"|TE.Density.Means$TE_type=="Unknown",]
#Subset to look only at TE families and rename column to TE_family
TE.Family.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="Unknown_fam"|TE.Density.Means$TE_type=="PIF_Harbinger"|TE.Density.Means$TE_type=="None"|TE.Density.Means$TE_type=="MULE"|TE.Density.Means$TE_type=="Copia"|TE.Density.Means$TE_type=="Gypsy"|TE.Density.Means$TE_type=="hAT"|TE.Density.Means$TE_type=="CMC_EnSpm",]
colnames(TE.Family.Density.Means)[2]<-"TE_Family"

#plot by TE type
ggplot(TE.Type.Density.Means)+geom_point(aes(x=Distance,y=avg,colour=TE_type,size=5))+scale_colour_manual(values=myColors)+facet_grid(~Location)

#plot by TE Family
ggplot(TE.Family.Density.Means)+geom_point(aes(x=Distance,y=avg,colour=TE_Family,size=5))+scale_colour_manual(values=myColors)+facet_grid(~Location)

write.csv(Reshaped,'data.csv')


