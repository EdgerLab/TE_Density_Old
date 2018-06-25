# By Scott Teresi
# The below script calculates TPM for a strawberry
# The code first has to create the data structure for the genes
#==================================================================
import sys, re, os, time, csv
from multiprocessing import Process
from collections import deque
#sys.path.insert(0, '/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/') # specify path for import
#import genic_elements
genes = deque() # use a deque structure for containing my elements
count_list = ['ERR855501.htseq.count','ERR855502.htseq.count','ERR855503.htseq.count','ERR855504.htseq.count','ERR855505.htseq.count','ERR855506.htseq.count']
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

def count_iterator(countfiles):
    a_dict = {}
    for countfile in countfiles:
        with open(countfile, 'r') as f_in:
            total_count = 0
            for row in f_in:
                row = re.split('\t+',row)
                key = row[0]
                count = int(row[1])

                if key not in a_dict:
                    a_dict[key] = count
                elif key in a_dict:
                    a_dict[key] += count
                total_count += count

    for countfile in countfiles:
        for elem in genes:
            name = elem.getName()
            count = a_dict[name]

            if countfile == countfiles[0]:
                elem.count_0 = count
            elif countfile == countfiles[1]:
                elem.count_1 = count
            elif countfile == countfiles[2]:
                elem.count_2 = count
            elif countfile == countfiles[3]:
                elem.count_3 = count
            elif countfile == countfiles[4]:
                elem.count_4 = count
            elif countfile == countfiles[5]:
                elem.count_5 = count

    TPM(total_count)
    #FPKM(total_count)

def TPM(total_count):
    for elem in genes:
        length = elem.length/1000 # we want it in kb
        for x in ['count_0','count_1','count_2','count_3','count_4','count_5']:
            count = getattr(elem,x)
            count_over_length = count/length

            #summation =
            setattr(elem,x,count)





    #for elem in genes:
        #for x in ['start', 'stop']:
            #print(getattr(elem, x))
            #print()



    #for elem in genes:
        #print(elem.__dict__)


        #b = [a for a in dir(elem) if not a.startswith('__')]
        #print(b)
        #for item in b:
            #if item[0:5] == 'count':
                #getattr(
                #item += 5







#=============================================================
def write_structure():
    with open(count_output, 'w') as f_out:
        fieldnames = ['name', 'chromosome', 'start', 'stop', 'length','count_0','count_1','count_2','count_3','count_4','count_5']
        f_out = csv.DictWriter(f_out,fieldnames=fieldnames)
        f_out.writeheader()
        for elem in genes:
            f_out.writerow(elem.__dict__)

def info():
    print('module name:', __name__)
    print('Started algorithm')
    print('parent process:', os.getppid())
    print('process id:', os.getpid())
    print()

def run_all(gtf_inputfile,mRNA_inputfile):
    info()
    gene_handler(mRNA_inputfile)
    gene_handler_2(gtf_inputfile)
    count_iterator(count_list)
    write_structure()

#=============================================================
class Gene(object):
    def __init__(self, maker_name, chromosome, start, stop ):
        self.name = maker_name
        self.chromosome = chromosome
        self.start = int(start)
        self.stop = int(stop)

        #count_list = ['ERR855501.htseq.count','ERR855502.htseq.count','ERR855503.htseq.count','ERR855504.htseq.count','ERR855505.htseq.count','ERR855506.htseq.count']
        self.count_0 = 0
        self.count_1 = 0
        self.count_2 = 0
        self.count_3 = 0
        self.count_4 = 0
        self.count_5 = 0

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
    #run_all('camarosa_gtf_data.gtf','mRNA.bed')
    run_all('xx00','mRNA00')



    #my_inputs = [['xx00','mRNA00']]#,['xx01','mRNA01'],['xx02','mRNA02'],['xx03','mRNA03'],['xx04','mRNA04'],['xx05','mRNA05'],
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


