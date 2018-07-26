# This separator file is different from the others because the order of the chromosomes is weird.
# Be careful when coding that you are working with the right sequence.

csplit H4_Genes.gtf /Fvb2/ /Fvb3/ /Fvb4/ /Fvb5/ /Fvb6/ /Fvb7/ -f gtf
echo ''
echo 'Starting mRNA'
echo ''
csplit mRNA.bed /Fvb2/ /Fvb3/ /Fvb4/ /Fvb5/ /Fvb6/ /Fvb7/ -f mRNA # -f mRNA creates a prefix of mRNA
