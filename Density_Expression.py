# By Scott Teresi
#----------------------------------------------------------

from Genic_Elements import Genic_Element, TE, Gene
import re
from collections import deque # I am going to implement a deque so that I can efficiently add TEs to my data set.
    # It won't be much use for lookup, not any better than a list, but it is good for constructing the set.
from Density_Algorithm import *
import csv
import time
import os
from multiprocessing import Process

#----------------------------------------------------------
transposons = deque()
genes = deque()
#----------------------------------------------------------
def te_handler(gff_inputfile):
    # One important caveat, I don't think I am grabbing the '1' class, of which there are 27 occurences, check with Pat.
    with open(gff_inputfile, 'r') as f_in:
        number = 1
        for row in f_in:
            row = re.split('\t+', row)
            if row[0][:10] == '##sequence':
                continue # remove weird sequence blocks
            classification = str(row[2])
            te_type = classification.split('/')[0]
            if te_type == 'unknown':
                te_type = 'Unknown'
            try:
                family = classification.split('/')[1]
            except IndexError:
                family = 'Unknown'
            if te_type == 'Simple_repeat':
                continue # removes the simple_repeats
            te_type = whitelist_type(te_type)
            family = whitelist_fam(family)

            if te_type == 'LTR' and family == 'Unknown':
                family = 'LTRUnknown'

            start = row[3]
            chromosome = str(row[0])
            stop = row[4]
            confidence = row[5]
            istart = int(start)
            istop = int(stop)
            length = istop - istart + 1
            name = str('TE_' + str(number))
            name = TE(name,chromosome,start,stop,length, te_type,family)
            transposons.append(name)
            number += 1


def whitelist_type(te_type):
    """Renames types into the correct types that we want, 'whitelisting' the acceptable types."""
    U = 'Unknown'
    N = 'None'
    master_type = {'RC?':'DNA','RC':'DNA', 'SINE?':U}
    for key, val in master_type.items():
        if te_type == key:
            te_type = val
    return te_type # returns te_type even if there is no match to anything

def whitelist_fam(family):
    """Renames families into the correct families that we want, 'whitelisting' the acceptable types."""
    U = 'Unknown'
    N = 'None'
    master_fam =  {'Uknown':U, 'MuDr':'MULE', 'MULE-MuDR':'MULE', 'Pao':U, 'Caulimovirus':U,
                   'hAT-Tag1':'hAT', 'hAT-Tip100':'hAT', 'hAT-Charlie':'hAT', 'Helitron':U, 'unknown':U,                                       'Maverick':U, 'Harbinger':'PIF-Harbinger', 'TcMar-Pogo':U, 'CR1':'LINE', 'hAT-Ac':'hAT',
                   'L2':'LINE', 'L1':'LINE', 'Jockey':'LINE', 'MuLE-MuDR':'MULE', 'MuDR':'MULE', 'Mutator':'MULE', 'Micro_like':U}
    for key,val in master_fam.items():
        if family == key:
            family = val
    return family

def gene_handler(mRNA_inputfile):
   """gene_handler() adds the genes into the 'genes' deque. Important to note that exon lengths come from elsewhere."""
   number = 1
   with open(mRNA_inputfile, 'r') as f_in:
      for row in f_in:
         row = re.split('\t+', row)
         chromosome = str(row[0])
         if chromosome[:6] == "contig":
            continue
         start = row[1]
         stop = row[2]
         maker_name = row[3].strip('\n')
         genes.append(Gene(number,chromosome,start,stop,maker_name))
         number +=1

def gene_handler_2(gtf_inputfile):
    """Gives us the exon lengths for each gene"""
    if selection == H4:
        a_dict = {}
        with open(gtf_inputfile,'r') as f_in:
            for row in f_in:
                row = re.split('\t+',row)
                exon_length = int(row[4]) - int(row[3]) + 1
                chromosome = str(row[0])
                classification = str(row[2])
                constructor = str(row[1])

                if chromosome[:6] == "contig" or chromosome[:3] != "Fvb":
                    continue

                if classification == 'exon' and constructor == 'maker':
                    gene = re.split('ID=|-mRNA',row[8])
                    key = gene[1].split(":")[0]

                    if key not in a_dict:
                        a_dict[key] = exon_length
                    elif key in a_dict:
                        a_dict[key] += exon_length

        for elem in genes:
            try:
                name = elem.getMaker_Name()
                length = int(a_dict[name])
                elem.length = length
            except KeyError:
                raise KeyError

    if selection == Camarosa:
        a_dict = {}
        with open(gtf_inputfile,'r') as f_in:
            for row in f_in:
                row = re.split('\t+',row)
                exon_length = int(row[4]) - int(row[3]) + 1
                chromosome = str(row[0])
                classification = str(row[2])
                constructor = str(row[1])

                if chromosome[:6] == "contig" or chromosome[:3] != "Fvb":
                    continue

                if classification == 'exon' and constructor == 'maker':
                    gene = re.split('ID=|-mRNA',row[8])
                    key = gene[1]

                    if key not in a_dict:
                        a_dict[key] = exon_length
                    elif key in a_dict:
                        a_dict[key] += exon_length

        for elem in genes:
            try:
                name = elem.getMaker_Name()
                length = int(a_dict[name])
                elem.length = length
            except KeyError:
                raise KeyError

def info():
    print('Started algorithm')
    print('parent process:', os.getppid())
    print('process id:', os.getpid())
    print()

def run_all(gff_inputfile, mRNA_inputfile, gtf_inputfile, file_list):
    te_handler(gff_inputfile)
    gene_handler(mRNA_inputfile)
    gene_handler_2(gtf_inputfile)

    info()
    #stats()
    get_densities(genes,transposons,500,500,10000,gtf_inputfile,file_list)

#---------------------------------------------------------
def stats():
    LTR_Count = 0
    Unknown_LTR_Count = 0
    All_Unknown_fam = 0
    All_Unknown = 0
    for elem in transposons:
        if elem.getTe_Type() == 'LTR':
            LTR_Count += 1
            if elem.getFamily() == 'Unknown':
                Unknown_LTR_Count += 1
        if elem.getFamily() == 'Unknown':
            All_Unknown_fam += 1
        if elem.getTe_Type() == 'Unknown':
            All_Unknown += 1

    print('LTR_Count ' + str(LTR_Count))
    print('Unknown_LTR_Count ' + str(Unknown_LTR_Count))
    print('All_Unknown_fam ' + str(All_Unknown_fam))
    print('All_Unknown ' + str(All_Unknown))
#---------------------------------------------------------




if __name__ == '__main__':
    H4 = 1
    Camarosa = 2

    selection = Camarosa ### THIS IS THE LINE TO MAKE YOUR SELECTION

    if selection == H4:
        print("Running H4")
        my_inputs = [['gtf00','mRNA00'],['gtf01','mRNA01'],['gtf02','mRNA02'],['gtf03','mRNA03'],['gtf04','mRNA04'],['gtf05','mRNA05'],['gtf06','mRNA06']]
        gff_inputfile = 'H4_TEs.gff'

    if selection == Camarosa:
        print("Running Camarosa")
        my_inputs = [['xx00','mRNA00'],['xx01','mRNA01'],['xx02','mRNA02'],['xx03','mRNA03'],['xx04','mRNA04'],['xx05','mRNA05'],['xx06','mRNA06']]
        gff_inputfile = 'camarosa_gff_data.gff'

    file_list = []
    for item in my_inputs:
        file_list.append(item[0])

    p_list = []
    for f_name in my_inputs:
        p = Process(target=run_all, args=(gff_inputfile, f_name[1], f_name[0], file_list,))
        p.start()
        p_list.append(p)

    for p in p_list:
        current_pid = p.pid
        p.join(None)
        print("pid '{}' joined!".format(current_pid))
