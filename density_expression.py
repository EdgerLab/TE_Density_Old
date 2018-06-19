# By Scott Teresi
#----------------------------------------------------------
from genic_elements import Genic_Element, TE, Gene
import re
from collections import deque # I am going to implement a deque so that I can efficiently add TEs to my data set.
    # It won't be much use for lookup, not any better than a list, but it is good for constructing the set.
from density_algorithm import *
import csv
import time
import os
#import multiprocessing
from multiprocessing import Process
#import multiprocessing

#----------------------------------------------------------
gff_inputfile = 'camarosa_gff_data.gff'
#gtf_inputfile = 'camarosa_gtf_data.gtf'
#gene_inputfile = 'mRNA.bed'
#genes_renamed = 'Fxa_chr_rename.txt'
transposons = deque()
genes = deque()
#----------------------------------------------------------
def te_handler():
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
    return te_type

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
    # This will have every gene because it is coming from the mRNA.bed file which has everything.
    # Oh so the problem is that some of the gene elements don't have a length. Because there is a mismatch
    # mismatch between the mRNA and the gtf partitions

    number = 1
    with open(mRNA_inputfile, 'r') as f_in:
        for row in f_in:
            row = re.split('\t+', row)
            chromosome = str(row[0])
            start = row[1]
            stop = row[2]
            maker_name = row[3].strip('\n')
            genes.append(Gene(number,chromosome,start,stop,maker_name))
            number +=1

def gene_handler_2(gtf_inputfile):
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
            name = elem.getMaker_Name()
            length = int(a_dict[name])
            elem.length = length
        except KeyError:
            continue

def info():
    #print('module name:', __name__)
    print('Started algorithm')
    print('parent process:', os.getppid())
    print('process id:', os.getpid())
    print()

def run_all(gtf_inputfile,mRNA_inputfile):
    te_handler()
    gene_handler(mRNA_inputfile)
    gene_handler_2(gtf_inputfile)

    #if gtf_inputfile == 'xx02':
        #for elem in genes:
            #print(elem.__dict__)
            #print(elem.length)
    info()
    get_densities(genes,transposons,500,500,10000,gtf_inputfile)

#---------------------------------------------------------
def reader(my_deque):
    pass

if __name__ == '__main__':
    #p = Pool()
    #p.map(run_all,['gtf_chunks_aa.gtf','gtf_chunks_ab.gtf','gtf_chunks_ac.gtf','gtf_chunks_ad.gtf','gtf_chunks_ae.gtf','gtf_chunks_af.gtf'])
    #for item in my_inputs:
        #p = Process(target=run_all, args=(item,))
        #p.start()
    #for item in my_inputs:
        #gtf_inputfile = item
        #te_handler()
        #gene_handler()
        #gene_handler_2(gtf_inputfile)

        #p = Process(target=get_densities, args=(genes,transposons,500,500,10000,gtf_inputfile,))
        #p = Process(target=run_all,args=(gtf_inputfile,))
        #p.start()

    #my_inputs = ['gtf_chunks_aa.gtf','gtf_chunks_ab.gtf','gtf_chunks_ac.gtf','gtf_chunks_ad.gtf','gtf_chunks_ae.gtf','gtf_chunks_af.gtf']
    #my_inputs = ['xx00','xx01','xx02','xx03','xx04','xx05','xx06']
    # for some reason it isn't working on anything past the xx00
    #my_inputs = ['xx01','xx02','xx03','xx04','xx05','xx06']
    my_inputs = [['xx00','mRNA00'],['xx01','mRNA01'],['xx02','mRNA02'],['xx03','mRNA03'],['xx04','mRNA04'],['xx05','mRNA05'],
                ['xx06','mRNA06']]
    p_list = []
    for f_name in my_inputs:
        p = Process(target=run_all, args=(f_name[0],f_name[1],))
        p.start()
        p_list.append(p)

    for p in p_list:
        current_pid = p.pid
        p.join(None)
        print("pid '{}' joined!".format(current_pid))
