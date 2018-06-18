# README file for the Camarosa Strawberry TE Density Project
# By Scott Teresi
---

# Considerations:

If you are adapting this code for your own use there are several things you will need to change in order to make it work for your system, these are in no particular order:

1. The TE types and families created in **whitelist\_type** and **whitelist\_fam** are unique to my analyses. You will have to change the dictionary there to sort TE classifications a different way.

2. I split up my genes into several files grouped by their chromosomes. To avoid working with a difficult gene file format, I make use of an mRNA bedfile to make getting all of the gene information easier. Both the bedfile and the gtf/gff file for the genes need to be split along chromosome identities. I make use of the script **separator.sh** for this.

3. There can be cases where a TE completely encompasses a gene from left to right, this is most likely an artifact of the TE annotation software, to deal with this, I do not add the TE values to the current density, and I output it to an overlap file so we can look at it later. Unfortunately I have hardcoded in the naming scheme for the overlap file. If you want to edit this you need to make sure the list matches with the *gtf\_inputfile* argument in **get\_densities**.

4. **separator.sh** needs to be edited to have the correct order of chromosomes as they appear in your gene and mRNA.bed files. Just make sure the chromosome sequence is right.
