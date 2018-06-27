# By Scott Teresi
# The below script creates a count file in preparation for RPKM, FPKM, and TPM
#==================================================================
import sys, re, os, time, csv
from multiprocessing import Process
from collections import deque

genes = deque() # use a deque structure for containing my elements
count_dict = {'ERR855501.htseq.count':'count_0','ERR855502.htseq.count':'count_1','ERR855503.htseq.count':'count_2','ERR855504.htseq.count':'count_3','ERR855505.htseq.count':'count_4','ERR855506.htseq.count':'count_5'}

# Specify your files here!!! Especially your count file and it's associated attribute

gene_dictionary = {}
count_output = 'expression_out.csv'
#==================================================================

def gene_handler(mRNA_inputfile):
    """gene_handler() adds the genes into the 'genes' deque. Important to note that exon lengths come from elsewhere."""
    # This will have every gene because it is coming from the mRNA.bed file which has everything.
    with open(mRNA_inputfile, 'r') as f_in:
        for row in f_in:
            row = re.split('\t+', row)
            chromosome = str(row[0])
            start = row[1]
            stop = row[2]
            maker_name = row[3].strip('\n')
            genes.append(Gene(maker_name,chromosome,start,stop))
            #gene_dictionary[maker_name]=self
            #Gene(maker_name,chromosome,start,stop)

def gene_handler_2(gtf_inputfile):
    """Use a simple dictionary to get the exon lengths calculated."""
    # Takes a gtf file as an input and iterates over to aggregate the correct exon lengths.
    a_dict = {}
    with open(gtf_inputfile,'r') as f_in:
        for row in f_in:
            row = re.split('\t+',row)
            exon_length = int(row[4]) - int(row[3]) + 1
            classification = str(row[2])
            constructor = str(row[1])
            if classification == 'exon' and constructor == 'maker':
                gene = re.split('ID=|-mRNA',row[8])
                key = gene[1]

                if key not in a_dict:
                    a_dict[key] = exon_length
                elif key in a_dict:
                    a_dict[key] += exon_length
                else:
                    raise ValueError("I don't know")
    for elem in genes:
        try:
            name = elem.getName()
            length = int(a_dict[name])
            elem.length = length
        except KeyError:
            continue

def count_iterator(count_dict):
    """Adds count values into the dictionary"""
    for key, val in count_dict.items():
        with open(key, 'r') as f_in:
            for row in f_in:
                row = re.split('\t+',row)
                gene = row[0]
                count = int(row[1])
                setattr(gene_dictionary[gene],val,count) # use a dictionary to get reference to the object with no loop
                # val is a reference to the attribute
                # count is the number we set for that attribute

def total_count():
    """Calculates the total count and adjusts the attribute accordingly"""
    for key, val in count_dict.items():
        total = 0
        for elem in genes:
            total += getattr(elem,val) # add the count value to the total
        specific_attribute = 'total_'+val
        for elem in genes:
            setattr(elem,specific_attribute,total)

def TPM(total_count):
    for elem in genes:
        length = elem.length/1000 # we want it in kb
        for x in ['count_0','count_1','count_2','count_3','count_4','count_5']:
            count = getattr(elem,x)
            count_over_length = count/length
            setattr(elem,x,count)

#=============================================================
def write_structure(a_dictionary):
    with open(count_output, 'w') as f_out:
        fieldnames = ['name', 'chromosome', 'start', 'stop', 'length','count_0','count_1','count_2','count_3','count_4','count_5', 'total_count_0','total_count_1','total_count_2','total_count_3','total_count_4','total_count_5'] # make sure your fieldnames contains strings of all your attributes
        f_out = csv.DictWriter(f_out,fieldnames=fieldnames)
        f_out.writeheader()
        for elem in a_dictionary:
            f_out.writerow(elem.__dict__)

def run_all(gtf_inputfile,mRNA_inputfile):
    gene_handler(mRNA_inputfile)
    gene_handler_2(gtf_inputfile)
    count_iterator(count_dict)
    total_count()
    write_structure(genes)

#=============================================================
class Gene(object):
    def __init__(self, maker_name, chromosome, start, stop ):
        self.name = maker_name
        self.chromosome = chromosome
        self.start = int(start)
        self.stop = int(stop)
        gene_dictionary[maker_name] = self # add the object reference into the dictionary

        # Make sure you get the right number of count and total count attributes made. I had 0-5

        self.count_0 = 0
        self.count_1 = 0
        self.count_2 = 0
        self.count_3 = 0
        self.count_4 = 0
        self.count_5 = 0

        self.total_count_0 = 0
        self.total_count_1 = 0
        self.total_count_2 = 0
        self.total_count_3 = 0
        self.total_count_4 = 0
        self.total_count_5 = 0

    def getName(self):
        return self.name

    def getChromosome(self):
        return self.chromosome

    def getStart(self):
        return self.start

    def getStop(self):
        return self.stop

    def getLength(self):
        return self.length

#===============================================================
if __name__ == '__main__':
    run_all('camarosa_gtf_data.gtf','mRNA.bed')
