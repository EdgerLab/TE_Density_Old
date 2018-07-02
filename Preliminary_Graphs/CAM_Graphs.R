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
TE_Data <- list.files(path=directory, pattern="*density_data.csv") # set the path to the current wd, and then grab all files with the cleaned name
rm(directory)
TE_Data <- do.call(rbind,lapply(TE_Data,read.csv)) # sometimes this step needs to be run twice?
TE_Data <- select(TE_Data,-number) # remove the numbers column, it is vestigial
#==========================================================================================================
Fvb1 = TE_Data[grep("Fvb1-*",TE_Data$chromosome),]
Fvb2 = TE_Data[grep("Fvb2-*",TE_Data$chromosome),]
Fvb3 = TE_Data[grep("Fvb3-*",TE_Data$chromosome),]
Fvb4 = TE_Data[grep("Fvb4-*",TE_Data$chromosome),]
Fvb5 = TE_Data[grep("Fvb5-*",TE_Data$chromosome),]
Fvb6 = TE_Data[grep("Fvb6-*",TE_Data$chromosome),]
Fvb7 = TE_Data[grep("Fvb7-*",TE_Data$chromosome),]
Fvb8 = TE_Data[grep("Fvb8-*",TE_Data$chromosome),]

# Re-shape data to make column for variable TE type/family and corresponding value for TE density
Reshaped<-melt(data = TE_Data,id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length","window_size"))
rm(TE_Data) # remove TE_Data because we need the RAM









# calculate the average and also divide by 2 to simplify, 2 used to mean full density.
TE.Density.Means<-Reshaped %>% group_by(variable,window_size) %>% summarise(avg=mean(value/2)) %>% arrange(avg) 

rm(Reshaped) # remove Reshaped because we need the RAM

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

# make left window negative
TE.Density.Means[TE.Density.Means$Location=="left",]$window_size<-TE.Density.Means[TE.Density.Means$Location=="left",]$window_size*-1

# Rename the middle column to have a capital I
TE.Density.Means[TE.Density.Means$Location=="intra",]$window_size<-"0"

# Subset to look only at TE types
TE.Type.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="DNA"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="LTR"|TE.Density.Means$TE_type=="Unknown",]

#Subset to look only at TE families and rename column to TE_family
# Will need to be renamed in a spot once the LINE fam is added
TE.Family.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="Unknown_fam"|TE.Density.Means$TE_type=="PIF_Harbinger"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="MULE"|TE.Density.Means$TE_type=="Copia"|TE.Density.Means$TE_type=="Gypsy"|TE.Density.Means$TE_type=="hAT"|TE.Density.Means$TE_type=="CMC_EnSpm",]
colnames(TE.Family.Density.Means)[2]<-"TE_Family" # LINE fam not none

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
# Plot along chromosome lines
# TE type
# Subset to look only at TE types

Reshaped<-melt(data = Fvb1,id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length","window_size"))
#rm(TE_Data) # remove TE_Data because we need the RAM









# calculate the average and also divide by 2 to simplify, 2 used to mean full density.
TE.Density.Means<-Reshaped %>% group_by(variable,window_size) %>% summarise(avg=mean(value/2)) %>% arrange(avg) 

rm(Reshaped) # remove Reshaped because we need the RAM

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

# make left window negative
TE.Density.Means[TE.Density.Means$Location=="left",]$window_size<-TE.Density.Means[TE.Density.Means$Location=="left",]$window_size*-1

# Rename the middle column to have a capital I
TE.Density.Means[TE.Density.Means$Location=="intra",]$window_size<-"0"

# Subset to look only at TE types
TE.Type.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="DNA"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="LTR"|TE.Density.Means$TE_type=="Unknown",]

#Subset to look only at TE families and rename column to TE_family
# Will need to be renamed in a spot once the LINE fam is added
TE.Family.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="Unknown_fam"|TE.Density.Means$TE_type=="PIF_Harbinger"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="MULE"|TE.Density.Means$TE_type=="Copia"|TE.Density.Means$TE_type=="Gypsy"|TE.Density.Means$TE_type=="hAT"|TE.Density.Means$TE_type=="CMC_EnSpm",]
colnames(TE.Family.Density.Means)[2]<-"TE_Family" # LINE fam not none

ggplot(TE.Type.Density.Means,aes(x=as.numeric(window_size),y=avg,group=TE_type,colour=TE_type))+geom_point(aes(size=5)) +
geom_line()+
scale_colour_manual(values=myColors)+facet_grid(~Location,scales="free") +
ylim(0,0.075) + ylab('TE Density') + xlab('Window Size')+scale_x_continuous(breaks = c(-10000, -9500, -9000, -8500, -8000, -7500, -7000, -6500, -6000, -5500, -5000, -4500, -4000, -3500, -3000, -2500, -2000, -1500,-1000,-500,0,500,1000,1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500, 10000))





#==============================================
# Expression Data
#library("stringi")
#library("tidyr")
#library("dplyr")
#library("reshape2")
#library("cowplot")

#directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/CAMDATA/")
#all_data <- list.files(path=directory, pattern="*density_data.csv") # set the path to the current wd, and then grab all files with the cleaned name
#all_data_df <- do.call(rbind,lapply(all_data,read.csv))  #put all the data into one dataframe
#all_data_df <- subset(all_data_df, select = -c(number)) # removes the number column
#Reshaped<-melt(data = all_data_df,id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length","window_size")) # only the values not mentioned, are the ones that turn into the variable and value column
#TE.Density.Means<-Reshaped %>% group_by(variable,window_size) %>% summarise(avg=mean(value)/2) %>% arrange(avg)

#-------------------
# Proximity
  library(scales)
  y_limit <- 100000
  y_step = 10000

  left_prox_graph <- ggplot(all_data_df, aes(x=prox_left)) +geom_histogram(binwidth=100,color='black',fill='white') +
  scale_x_continuous(breaks=seq(0,5000,500),limits = c(0, 5000)) +
  scale_y_continuous(breaks = seq(0,y_limit,y_step),limits=c(0,y_limit),labels=comma) +
  labs(title='Left Side', x = 'BP Away to Closest TE', y = 'Number of Genes')
  print(left_prox_graph)

  right_prox_graph <- ggplot(all_data_df, aes(x=prox_right)) +geom_histogram(binwidth=100,color='black',fill='white') +
  scale_x_continuous(breaks=seq(0,5000,500),limits = c(0, 5000)) +
  scale_y_continuous(breaks = seq(0,y_limit,y_step),limits=c(0,y_limit),labels=comma) +
  labs(title='Right Side', x = 'BP Away to Closest TE', y = 'Number of Genes')
  print(right_prox_graph)
   #Aesthetics just tell what the x and y axis should be, coloring and such
   #good habit is to put as much general aesthetics into the main ggplot function



#============================================================================

directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/CAMDATA/")
TE_Data <- list.files(path=directory, pattern="*density_data.csv") # set the path to the current wd, and then grab all files with the cleaned name
TE_Data <- do.call(rbind,lapply(TE_Data,read.csv))
TE_Data <- select(TE_Data,-number) # remove the numbers column, it is vestigial

second_directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/Expression/") # specify the file path
Expression <- read.csv('expression_out.csv', header=TRUE, sep=",") # read the expression csv

Expression$TPM_AVG <- rowMeans(Expression[c('TPM_0','TPM_1','TPM_2','TPM_3','TPM_4','TPM_5')]) # calculate the mean for TPM
Expression$FPKM_AVG <- rowMeans(Expression[c('FPKM_0','FPKM_1','FPKM_2','FPKM_3','FPKM_4','FPKM_5')]) # calculate the mean for FPKM

GeneExpressionDF<-Expression %>% select(name,TPM_AVG,FPKM_AVG) # make a new data frame with only the TPM and FPKM averages
colnames(GeneExpressionDF)<-c("maker_name","TPM_AVG","FPKM_AVG") # reset the name column to be maker_name
Reshapetest<-Reshaped %>%  left_join(GeneExpressionDF, by="maker_name")

TE.TypeExp<-Reshapetest %>% select(one_of(c("DNA_left","DNA_right","LINE_left","LINE_right","LTR_left","LTR_right","Unknown_left","Unknown_right")))
#Subset to look only at TE families and rename column to TE_family

# Will need to be renamed in a spot once the LINE fam is added
#TE.FamilyExp<-TE.Density.Means[TE.Density.Means$TE_type=="Unknown_fam"|TE.Density.Means$TE_type=="PIF_Harbinger"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="MULE"|TE.Density.Means$TE_type=="Copia"|TE.Density.Means$TE_type=="Gypsy"|TE.Density.Means$TE_type=="hAT"|TE.Density.Means$TE_type=="CMC_EnSpm",]


#ggsave(ggplot(Reshapetest,aes(x=window_size,y=TPM_AVG),groupby=variable,colour=variable)+geom_point()+facet_wrap(~variable))


#==============================================
# Write the data
#write.csv(Reshaped,'data.csv')

