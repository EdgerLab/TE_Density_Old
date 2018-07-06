library("stringi")
library("tidyr")
library("dplyr")
library("reshape2")
library("cowplot")
directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/CAMDATA/")
TE_Data <- list.files(path=directory, pattern="*density_data.csv") # set the path to the current wd, and then grab all files with the cleaned name
TE_Data <- do.call(rbind,lapply(TE_Data,read.csv))
TE_Data <- select(TE_Data,-number) # remove the numbers column, it is vestigial

#==========================================================================================================


# Re-shape data to make column for variable TE type/family and corresponding value for TE density
TE_Data = melt(data = TE_Data,id=c("chromosome","maker_name","start","stop","prox_left","prox_right","they_are_inside","length","window_size"))
TE_Data$value = TE_Data$value/2
TE_Data = TE_Data[TE_Data$value<=1,]

directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/GO/")
write.csv(TE_Data,file = 'TE_Data.csv')

group_by(TE_Data,window_size) %>% do(data.frame(t(quantile(.$value, probs =  c(0.1,0.99)))))
window_groupings = group_by(TE_Data,window_size,variable) %>% do(data.frame(t(quantile(.$value, probs =  c(0.1,0.9)))))

directory = setwd("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/GO/")
write.csv(window_groupings,file = 'window_groupings.csv')

#TE_Data_2 = tbl_df(TE_Data)
#output = TE_Data[TE_Data$window_size==500,]#,TE_Data$variable=='DNA_left']
#output = TE_Data[TE_Data$window_size==500 & TE_Data$variable=='DNA_left' & TE_Data$value > 0.133,]#,TE_Data$variable=='DNA_left']

