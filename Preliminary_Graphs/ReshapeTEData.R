# By Scott Teresi and Kevin Bird
#==========================================================================================================
#install.packages("reshape2")
#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("cowplot")
#==========================================================================================================
library("stringi")
library("tidyr")
library("dplyr")
library("reshape2")
library("cowplot")
#==========================================================================================================
# Load the data
all_data <- list.files(path=getwd(), pattern="*cleaned_data.csv") # set the path to the current wd, and then grab all files with the cleaned name
all_data_df <- do.call(rbind,lapply(all_data,function(x){read.csv(x,quote="'")}))  #put all the data into one dataframe

#`500_cleaned_data` <- read.csv("500_xx00_cleaned_data.csv", quote="'")
#==========================================================================================================
# Re-shape data to make column for variable TE type/family and corresponding value for TE density
Reshaped<-melt(data = all_data_df[-1],id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length","windowSize"))
TE.Density.Means<-Reshaped %>% group_by(variable,windowSize) %>% summarise(avg=mean(value)/2) %>% arrange(avg)
# Split TE type/family and location (e.g. left, right, intragenic) in hacky way (reverse string, split by "_", reverse strings again)
TE.Density.Means$variable <-stri_reverse(TE.Density.Means$variable)
TE.Density.Means<-TE.Density.Means%>%separate(variable,into = c("Location","TE_type"),sep="_",extra="merge")
TE.Density.Means$Location<-stri_reverse(TE.Density.Means$Location)
TE.Density.Means$TE_type<-stri_reverse(TE.Density.Means$TE_type)

# add Column for window size
#TE.Density.Means$Distance<-"500bp" # this is a manual old step

# Make colourblind friendly palette :)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
myColors <- colorRampPalette(cbPalette)(12)
names(myColors) <- unique(TE.Density.Means$TE_type)
colScale <- scale_fill_manual(name="", values=myColors)

# Reorder Location factors so that plot shows left, intra, right
TE.Density.Means$Location <- factor(TE.Density.Means$Location, levels=c("left", "intra", "right"))

TE.Density.Means[TE.Density.Means$Location=="left",]$windowSize<-TE.Density.Means[TE.Density.Means$Location=="left",]$windowSize*-1

# Rename the middle column to have a capital I
TE.Density.Means[TE.Density.Means$Location=="intra",]$windowSize<-"0"
#TE.Density.Means$windowSize <- factor(levels=c(500,1000,1500))
#TE.Density.Means$windowSize<-as.factor(TE.Density.Means$windowSize)
#TE.Density.Means$windowSize <- factor(TE.Density.Means$windowSize, levels=c(-1500,-1000,-500,500, 1000, 1500,"Intra_Gene"))
#TE.Density.Means[TE.Density.Means$Location=="intra"$windowSize]$levels=c(1500,1000,500)

# Make columns
#TE.Density.Means$Line_Values <- '' # made a empty line
#TE.Density.Means[TE.Density.Means$windowSize==500,]$Line_Values <- 1
#TE.Density.Means[TE.Density.Means$windowSize==1000,]$Line_Values <- 2
#TE.Density.Means[TE.Density.Means$windowSize==1500,]$Line_Values <- 3 
#TE.Density.Means[TE.Density.Means$windowSize=='Intra_Gene',]$Line_Values <- 0 

#TE.Density.Means[TE.Density.Means$windowSize==-500,]$Line_Values <- 1
#TE.Density.Means[TE.Density.Means$windowSize==-1000,]$Line_Values <- 2
#TE.Density.Means[TE.Density.Means$windowSize==-1500,]$Line_Values <- 3 
#TE.Density.Means[TE.Density.Means$windowSize=='Intra_Gene',]$Line_Values <- 0 

# Subset to look only at TE types
TE.Type.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="DNA"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="LTR"|TE.Density.Means$TE_type=="Unknown",]

#Subset to look only at TE families and rename column to TE_family
# !!!!!!!!!!!!

# Will need to be renamed in a spot once the LINE fam is added
TE.Family.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="Unknown_fam"|TE.Density.Means$TE_type=="PIF_Harbinger"|TE.Density.Means$TE_type=="None"|TE.Density.Means$TE_type=="MULE"|TE.Density.Means$TE_type=="Copia"|TE.Density.Means$TE_type=="Gypsy"|TE.Density.Means$TE_type=="hAT"|TE.Density.Means$TE_type=="CMC_EnSpm",]
colnames(TE.Family.Density.Means)[2]<-"TE_Family"


#==============================================
# Plot
# TE type
ggplot(TE.Type.Density.Means,aes(x=as.numeric(windowSize),y=avg,group=TE_type,colour=TE_type))+geom_point(aes(size=5)) + 
geom_line()+
scale_colour_manual(values=myColors)+facet_grid(~Location,scales="free") + 
ylim(0,1) + ylab('TE Density') + xlab('Window Size')+scale_x_continuous(breaks = c(-1500,-1000,-500,0,500,1000,1500))

# TE Family
ggplot(TE.Family.Density.Means,aes(x=as.numeric(windowSize),y=avg,group=TE_Family,colour=TE_Family))+geom_point(aes(size=5)) + 
geom_line()+
scale_colour_manual(values=myColors)+facet_grid(~Location,scales="free") + 
ylim(0,1) + ylab('TE Density') + xlab('Window Size')+scale_x_continuous(breaks = c(-1500,-1000,-500,0,500,1000,1500))
#==============================================






#==============================================
# Write the data
#write.csv(Reshaped,'data.csv')