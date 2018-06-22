# By Scott Teresi
# The below script calculates TPM for a strawberry
#====================================
import sys, re, os, time, csv
from multiprocessing import Process
from collections import deque
sys.path.insert(0, '/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/') # specify path for import
import genic_elements
genes = deque()
#====================================

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
    print('module name:', __name__)
    print('Started algorithm')
    print('parent process:', os.getppid())
    print('process id:', os.getpid())
    print()

def run_all(gtf_inputfile,mRNA_inputfile):
    gene_handler(mRNA_inputfile)
    gene_handler_2(gtf_inputfile)

    #if gtf_inputfile == 'xx02':
        #for elem in genes:
            #print(elem.__dict__)
            #print(elem.length)
    info()


#===============================================================
if __name__ == '__main__':
    #run_all('xx00','mRNA00')
    print('Yeah')



    #my_inputs = [['xx00','mRNA00'],['xx01','mRNA01'],['xx02','mRNA02'],['xx03','mRNA03'],['xx04','mRNA04'],['xx05','mRNA05'],
                #['xx06','mRNA06']]
    #p_list = []
    #for f_name in my_inputs:
        #p = Process(target=run_all, args=(f_name[0],f_name[1],))
        #p.start()
        #p_list.append(p)

    #for p in p_list:
        #current_pid = p.pid
        #p.join(None)
        #print("pid '{}' joined!".format(current_pid))
