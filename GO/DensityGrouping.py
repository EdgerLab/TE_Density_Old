# By Scott Teresi --- July 2018
#---------------------------------------------
import csv
from collections import deque
import subprocess
#import glob

# This function groups the TE's of my file into bins for GO enrichment.
# I use this file to partition the density of my genes

#file_list = glob.glob("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/CAMDATA/*density_data.csv")
#---------------------------------------------
left = deque()
intra = deque()
right = deque()

def write_structure():
    with open('GO_Groupings.csv','w') as f_out:
        my_header = ['chromosome','name','start','stop','window_size','classification','value']
        writer = csv.writer(f_out)
        writer.writerow([i for i in my_header])

        with open("TE_Data.csv") as csv_file:
            reader = csv.DictReader(csv_file) # if fieldnames is not supplied to DictReader then it uses the first line to make the fieldnames
            for row in reader: # refer to the variables through string name
                name = row['maker_name']
                window_size = int(row['window_size'])
                classification = row['variable']

                if classification[:4] == 'LTR_': # don't want LTRUnknown

                    value = float(row['value'])
                    chromosome = row['chromosome']
                    start = row['start']
                    stop = row['stop']

                    combined_name = str(window_size)+ '_' + classification

                    #if value <= Lower[combined_name][2] or value >= Upper[combined_name][2]:
                        #write_this = [chromosome,name,start,stop,window_size,classification,value]
                        #writer.writerow(write_this)

                    if window_size == 1000 and classification[:3] == 'LTR':
                        if value >= float(Upper[combined_name][2]): # grab the upper bound values
                            element = Gene(chromosome,name,start,stop,window_size,classification,value)
                            #print(element.__dict__)

                            if classification[-4:] == 'left':
                                left.append(element)

                            elif classification[-5:] == 'right':
                                right.append(element)

                            elif classification[-5:] == 'intra':
                                intra.append(element)


    fieldnames = ['chromosome','name','start','stop','window_size','classification','value']
    with open('left_bounded.txt','w') as f_out:
        f_out = csv.DictWriter(f_out,fieldnames=fieldnames)
        f_out.writeheader()
        for elem in left:
            f_out.writerow(elem.__dict__)

    #bashcommand = 'sort -r -k7,7 -n -t"," left_bounded.txt > sort_left.txt'
    #process = subprocess.Popen(bashcommand.split(), stdout = subprocess.PIPE)
    #output,error = process.communicate()

    with open('intra_bounded.txt','w') as f_out:
        f_out = csv.DictWriter(f_out,fieldnames=fieldnames)
        f_out.writeheader()
        for elem in intra:
            f_out.writerow(elem.__dict__)
    #bashcommand = 'sort -r -k7,7 -n -t"," intra_bounded.txt > sort_intra.txt'
    #process = subprocess.Popen(bashcommand.split(), stdout = subprocess.PIPE)
    #output,error = process.communicate()

    with open('right_bounded.txt','w') as f_out:
        f_out = csv.DictWriter(f_out,fieldnames=fieldnames)
        f_out.writeheader()
        for elem in right:
            f_out.writerow(elem.__dict__)
    #bashcommand = 'sort -r -k7,7 -n -t"," right_bounded.txt > sort_right.txt'
    #process = subprocess.Popen(bashcommand.split(), stdout = subprocess.PIPE)
    #output,error = process.communicate()




def whitelist_maker():
    """Make a dictionary of the correct pairs"""
    global Lower; Lower = {}
    global Upper; Upper = {}
    with open('Percentile_Whitelist.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            window_size = int(row['window_size'])
            classification = str(row['variable'])
            name = str(window_size)+ '_' + classification
            lower_number = float(row['X10.'])
            upper_number = float(row['X90.'])

            Lower[name] = (window_size,classification,lower_number)
            Upper[name] = (window_size,classification,upper_number)

class Gene(object):
    def __init__(self, chromosome, maker_name, start, stop, window_size, classification, value ):
        self.name = maker_name
        self.chromosome = chromosome
        self.start = int(start)
        self.stop = int(stop)
        self.classification = classification
        self.window_size = window_size
        self.value = value
        #gene_dictionary[maker_name]=self



if __name__ == '__main__':
    #classification = 'DNA_left'
    #print(classification[-5:])
    whitelist_maker()
    write_structure()

