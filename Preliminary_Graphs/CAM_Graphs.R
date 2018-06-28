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
directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/CAMDATA/")
all_data <- list.files(path=directory, pattern="*density_data.csv") # set the path to the current wd, and then grab all files with the cleaned name
all_data_df <- do.call(rbind,lapply(all_data,read.csv))  #put all the data into one dataframe
all_data_df <- subset(all_data_df, select = -c(number))
#==========================================================================================================
# Re-shape data to make column for variable TE type/family and corresponding value for TE density

Reshaped<-melt(data = all_data_df,id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length","window_size"))
TE.Density.Means<-Reshaped %>% group_by(variable,window_size) %>% summarise(avg=mean(value)/2) %>% arrange(avg)
# Split TE type/family and location (e.g. left, right, intragenic) in hacky way (reverse string, split by "_", reverse strings again)
TE.Density.Means$variable <-stri_reverse(TE.Density.Means$variable)
TE.Density.Means<-TE.Density.Means%>%separate(variable,into = c("Location","TE_type"),sep="_",extra="merge")
TE.Density.Means$Location<-stri_reverse(TE.Density.Means$Location)
TE.Density.Means$TE_type<-stri_reverse(TE.Density.Means$TE_type)

# Make colourblind friendly palette :)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
myColors <- colorRampPalette(cbPalette)(12)
names(myColors) <- unique(TE.Density.Means$TE_type)
colScale <- scale_fill_manual(name="", values=myColors)

# Reorder Location factors so that plot shows left, intra, right
TE.Density.Means$Location <- factor(TE.Density.Means$Location, levels=c("left", "intra", "right"))

TE.Density.Means[TE.Density.Means$Location=="left",]$window_size<-TE.Density.Means[TE.Density.Means$Location=="left",]$window_size*-1

# Rename the middle column to have a capital I
TE.Density.Means[TE.Density.Means$Location=="intra",]$window_size<-"0"

# Subset to look only at TE types
TE.Type.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="DNA"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="LTR"|TE.Density.Means$TE_type=="Unknown",]

#Subset to look only at TE families and rename column to TE_family

# Will need to be renamed in a spot once the LINE fam is added
TE.Family.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="Unknown_fam"|TE.Density.Means$TE_type=="PIF_Harbinger"|TE.Density.Means$TE_type=="None"|TE.Density.Means$TE_type=="MULE"|TE.Density.Means$TE_type=="Copia"|TE.Density.Means$TE_type=="Gypsy"|TE.Density.Means$TE_type=="hAT"|TE.Density.Means$TE_type=="CMC_EnSpm",]
colnames(TE.Family.Density.Means)[2]<-"TE_Family"

#==============================================
# Plot
# TE type
ggplot(TE.Type.Density.Means,aes(x=as.numeric(window_size),y=avg,group=TE_type,colour=TE_type))+geom_point(aes(size=5)) +
geom_line()+
scale_colour_manual(values=myColors)+facet_grid(~Location,scales="free") +
ylim(0,0.075) + ylab('TE Density') + xlab('Window Size')+scale_x_continuous(breaks = c(-10000, -9500, -9000, -8500, -8000, -7500, -7000, -6500, -6000, -5500, -5000, -4500, -4000, -3500, -3000, -2500, -2000, -1500,-1000,-500,0,500,1000,1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000))

# TE Family
ggplot(TE.Family.Density.Means,aes(x=as.numeric(window_size),y=avg,group=TE_Family,colour=TE_Family))+geom_point(aes(size=5)) +
geom_line()+
scale_colour_manual(values=myColors)+facet_grid(~Location,scales="free") +
ylim(0,0.075) + ylab('TE Density') + xlab('Window Size')+scale_x_continuous(breaks = c(-10000, -9500, -9000, -8500, -8000, -7500, -7000, -6500, -6000, -5500, -5000, -4500, -4000, -3500, -3000, -2500, -2000, -1500,-1000,-500,0,500,1000,1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000))

#==============================================

#==============================================
# Write the data
#write.csv(Reshaped,'data.csv')
