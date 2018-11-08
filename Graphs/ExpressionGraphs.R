# By Scott Teresi
#------------------------------
library("stringi")
library("tidyr")
library("dplyr")
library("reshape2")
library("cowplot")
#==========================================================================================================
# Load the data
directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/CAMDATA/")
TE_Data <- list.files(path=directory, pattern="*density_data.csv") # set the path to the current wd, and then grab all files with the cleaned name
TE_Data <- do.call(rbind,lapply(TE_Data,read.csv))

TE_Data$DNAAVG <- TE_Data$DNA_left+TE_Data$DNA_right # do not divide by 2 yet
TE_Data$LTRAVG <- TE_Data$LTR_left+TE_Data$LTR_right # do not divide by 2 yet
#write.csv(TE_Data, file='TE_Data.csv')

#TE_Data <- select(TE_Data, -DNA_left)
#TE_Data <- select(TE_Data, -DNA_right)
#TE_Data <- select(TE_Data, -LTR_left)
#TE_Data <- select(TE_Data, -LTR_right)

#TE_Data[is.nan(TE_Data$LTRAVG),]$LTRAVG <- 0
#TE_Data[is.infinite(TE_Data$LTRAVG),]$LTRAVG <- 0
#TE_Data[is.nan(TE_Data$DNAAVG),]$DNAAVG <- 0
#TE_Data[is.infinite(TE_Data$DNAAVG),]$DNAAVG <- 0

TE_Data <- select(TE_Data, -number) # remove the numbers column, it is vestigial
colnames(TE_Data)[colnames(TE_Data)=="Unknown_fam_left"] <- "UnknownFam_left"
colnames(TE_Data)[colnames(TE_Data)=="Unknown_fam_intra"] <- "UnknownFam_intra"
colnames(TE_Data)[colnames(TE_Data)=="Unknown_fam_right"] <- "UnknownFam_right"

Reshaped<-melt(data = TE_Data,id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length","window_size"))
rm(TE_Data) # remove TE_Data because we need the RAM

# Remove the 0s from the data frame
Reshaped <- Reshaped[Reshaped$value!=0,]

#Reshaped$value <- Reshaped$value/2 # stop dividng my two
Reshaped = Reshaped[Reshaped$value<=1,]

#TE.Density.Means<-Reshaped %>% group_by(variable,window_size) %>% summarise(avg=mean(value)/2) %>% arrange(avg)
#Reshaped2 <- Reshaped %>% group_by(variable,window_size) %>% summarise(avg=mean(value)/2) %>% arrange(avg)



# Split TE type/family and location (e.g. left, right, intragenic) in hacky way (reverse string, split by "_", reverse strings again)
#TE.Density.Means$variable <-stri_reverse(TE.Density.Means$variable)
#TE.Density.Means<-TE.Density.Means%>%separate(variable,into = c("Location","TE_type"),sep="_",extra="merge")
#TE.Density.Means$Location<-stri_reverse(TE.Density.Means$Location)
#TE.Density.Means$TE_type<-stri_reverse(TE.Density.Means$TE_type)

# Make colourblind friendly palette :)
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#myColors <- colorRampPalette(cbPalette)(14)
#names(myColors) <- unique(TE.Density.Means$TE_type)
#colScale <- scale_fill_manual(name="", values=myColors)

# Reorder Location factors so that plot shows left, intra, right
#TE.Density.Means$Location <- factor(TE.Density.Means$Location, levels=c("left", "intra", "right"))

# Make left side windows negative
#TE.Density.Means[TE.Density.Means$Location=="left",]$window_size<-TE.Density.Means[TE.Density.Means$Location=="left",]$window_size*-1

# Rename the middle column to have a capital I
#TE.Density.Means[TE.Density.Means$Location=="intra",]$window_size<-"0"

# Subset to look only at TE types
#TE.Type.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="DNA"|TE.Density.Means$TE_type=="LINE"|TE.Density.Means$TE_type=="LTR"|TE.Density.Means$TE_type=="Unknown"|TE.Density.Means$TE_type=="LTRUnknown",]

#Subset to look only at TE families and rename column to TE_family

#TE.Family.Density.Means<-TE.Density.Means[TE.Density.Means$TE_type=="Unknown_fam"|TE.Density.Means$TE_type=="PIF_Harbinger"|TE.Density.Means$TE_type=="LINE_fam"|TE.Density.Means$TE_type=="MULE"|TE.Density.Means$TE_type=="Copia"|TE.Density.Means$TE_type=="Gypsy"|TE.Density.Means$TE_type=="hAT"|TE.Density.Means$TE_type=="CMC_EnSpm",]
#colnames(TE.Family.Density.Means)[2]<-"TE_Family" # LINE fam not none

#============================================================================
setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/Expression/test_graphs/") # specify the file path


Expression <- read.csv('expression_out.csv', header=TRUE, sep=",") # read the expression csv
Expression$TPM_AVG <- rowMeans(Expression[c('TPM_0','TPM_1','TPM_2','TPM_3','TPM_4','TPM_5')]) # calculate the mean for TPM
Expression$FPKM_AVG <- rowMeans(Expression[c('FPKM_0','FPKM_1','FPKM_2','FPKM_3','FPKM_4','FPKM_5')]) # calculate the mean for FPKM
GeneExpressionDF<-Expression %>% select(name,TPM_AVG,FPKM_AVG) # make a new data frame with only the TPM and FPKM averages
colnames(GeneExpressionDF)<-c("maker_name","TPM_AVG","FPKM_AVG") # reset the name column to be maker_name


#setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/Expression/test_graphs/") # specify the file path

# Take a subset of Reshaped by a variable
TE.TypeExp <-Reshaped %>% filter(variable == "LTR_left")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('LTR_left.pdf')
ggplot(TE.TypeExp,aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="LTR Tranposon Density Upstream", x = "Tranposon Density by Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()

# Take a subset of Reshaped by a variable
TE.TypeExp <-Reshaped %>% filter(variable == "LTR_right")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('LTR_right.pdf')
ggplot(TE.TypeExp,aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="LTR Tranposon Density Downstream", x = "Tranposon Density by Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()


# Take a subset of Reshaped by a variable
TE.TypeExp <-Reshaped %>% filter(variable == "LTRAVG")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('LTRAVG.pdf')
ggplot(TE.TypeExp,aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="Average LTR Transposon Density (Upstream and Downstream)", x = "Tranposon Density by Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()

#---------------------------------------------------------
# DNA


TE.TypeExp <-Reshaped %>% filter(variable == "DNA_left")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('DNA_left.pdf')
ggplot(TE.TypeExp,aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="DNA Tranposon Density Upstream", x = "Tranposon Density by Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()

# Take a subset of Reshaped by a variable
TE.TypeExp <-Reshaped %>% filter(variable == "DNA_right")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('DNA_right.pdf')
ggplot(TE.TypeExp,aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="DNA Tranposon Density Downstream", x = "Tranposon Density by Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()


# Take a subset of Reshaped by a variable
TE.TypeExp <-Reshaped %>% filter(variable == "DNAAVG")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('DNAAVG.pdf')
ggplot(TE.TypeExp,aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="Average DNA Transposon Density (Upstream and Downstream)", x = "Tranposon Density by Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()

#--------------------------------------------
# LTR sub-windows

TE.TypeExp <-Reshaped %>% filter(variable == "LTR_left")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
#TE.TypeExp <- TE.TypeExp[TE.TypeExp$TPM_AVG>=0.1,] # filter out things with an expression < 0.1
pdf('LTR_left_window.pdf')
ggplot(TE.TypeExp[TE.TypeExp$window_size==1000,],aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_point(position='jitter')+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="LTR Tranposon Density Upstream", x = "Tranposon Density in a Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
#scale_y_continuous(limits=c(-39,-43))
dev.off()

# Take a subset of Reshaped by a variable
TE.TypeExp <-Reshaped %>% filter(variable == "LTR_right")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
#TE.TypeExp <- TE.TypeExp[TE.TypeExp$TPM_AVG>=0.1,] # filter out things with an expression < 0.1
pdf('LTR_right_window.pdf')
ggplot(TE.TypeExp[TE.TypeExp$window_size==1000,],aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_point(position='jitter')+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="LTR Transposon Density Downstream", x = "Tranposon Density in a Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
#scale_y_continuous(limits=c(-39,-43))
dev.off()

TE.TypeExp <-Reshaped %>% filter(variable == "LTRAVG")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
#TE.TypeExp <- TE.TypeExp[TE.TypeExp$TPM_AVG>=0.1,] # filter out things with an expression < 0.1
pdf('LTRAVG_window.pdf')
ggplot(TE.TypeExp[TE.TypeExp$window_size==1000,],aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_point(position='jitter')+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="Average LTR Transposon Density (Upstream and Downstream)", x = "Tranposon Density in a Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
#ylim(-39,-43)
#scale_y_continuous(limits=(-39,-43))
dev.off()

#---------------------------------------------
# DNA sub-windows

TE.TypeExp <-Reshaped %>% filter(variable == "DNA_left")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('DNA_left_window.pdf')
ggplot(TE.TypeExp[TE.TypeExp$window_size==1000,],aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="DNA Transposon Density Upstream", x = "Tranposon Density in a Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()

# Take a subset of Reshaped by a variable
TE.TypeExp <-Reshaped %>% filter(variable == "DNA_right")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('DNA_right_window.pdf')
ggplot(TE.TypeExp[TE.TypeExp$window_size==1000,],aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="DNA Transposon Density Downstream", x = "Tranposon Density in a Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()

TE.TypeExp <-Reshaped %>% filter(variable == "DNAAVG")
TE.TypeExp <- TE.TypeExp %>%  left_join(GeneExpressionDF, by="maker_name")
TE.TypeExp <- TE.TypeExp[TE.TypeExp$value<=1,] # filter out things with a density above 1
pdf('DNAAVG_window.pdf')
ggplot(TE.TypeExp[TE.TypeExp$window_size==1000,],aes(x=value,y=log2(TPM_AVG),color=as.character(window_size)))+geom_smooth(method='auto')+facet_wrap(~window_size)+
labs(title="Average DNA Transposon Density (Upstream and Downstream)", x = "Tranposon Density in a Window", y= "Gene Expression (Log(2) TPM)", color="Window Size")+
scale_x_continuous(breaks=seq(0,1))
dev.off()




#----------------------------------------------
# Put em together+
#plot_grid(DNARIGHT_window_plot,DNALEFT_window_plot,LTRLEFT_window_plot,LTRRIGHT_window_plot,nrow=2,scale=0.8)#labels=c('DNA Transposon Density Downstream','DNA Transposon Density Upstream', 'LTR Transposon Density Upstream', 'LTR Transposon Density Downstream'),nrow=2)



