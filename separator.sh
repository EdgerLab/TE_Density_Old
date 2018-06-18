#split -l 253729 --additional-suffix=.gtf camarosa_gtf_data.gtf gtf_chunks_
csplit camarosa_gtf_data.gtf /Fvb2/ /Fvb3/ /Fvb4/ /Fvb5/ /Fvb6/ /Fvb7/ #-f gtf
#csplit camarosa_gff_data.gff /Fvb2/ /Fvb3/ /Fvb4/ /Fvb5/ /Fvb6/ /Fvb7/ -f gff
csplit mRNA.bed /Fvb2/ /Fvb3/ /Fvb4/ /Fvb5/ /Fvb6/ /Fvb7/ -f mRNA
#------------------------------------------------------------
#lines=(`wc -l camarosa_gtf_data.gtf | cut -f1 -d' '`)
#echo $lines
#numbers_half=$(( lines / 2))
#echo $numbers_half
#numbers_quarter=$(( numbers_half / 2))
#numbers_quarter=$((numbers_quarter+1))
#echo $numbers_quarter
# I need to add one ???
#-------------------
#head camarosa_gtf_data.gtf -n $numbers_half > First_Half_Genes.gtf
#tail camarosa_gtf_data.gtf -n $numbers_half > Second_Half_Genes.gtf

#head First_Half_Genes.gtf -n $numbers_quarter > First_Quarter_Genes.gtf
#tail First_Half_Genes.gtf -n $numbers_quarter > Second_Quarter_Genes.gtf

#head Second_Half_Genes.gtf -n $numbers_quarter > Third_Quarter_Genes.gtf
#tail Second_Half_Genes.gtf -n $numbers_quarter > Fourth_Quarter_Genes.gtf

#rm -r {First_Half_Genes.gtf,Second_Half_Genes.gtf}
