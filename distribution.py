# By Scott Teresi
#----------------------------------------------------------
import re
from collections import deque # I am going to implement a deque so that I can efficiently add TEs to my data set.
    # It won't be much use for lookup, not any better than a list, but it is good for constructing the set.
from density_algorithm import *
import csv
import time
import os
from multiprocessing import Process
#----------------------------------------------------------
gff_inputfile = 'camarosa_gff_data.gff'
transposons = deque()
TE_Table = 'TE_Table.csv'
Type_Table = 'Type_Table.csv'
Family_Table = 'Family_Table.csv'
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
    distribution()


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


def distribution():
    with open(TE_Table, 'w') as f_out:
        fieldnames = ['number', 'chromosome', 'start', 'stop', 'length', 'te_type', 'family'] # Supply the header to the file, this also orders your output columns in the supplied order. Must have all attributes listed.
        # You will need to edit the fieldnames to match your data
        # You may not leave out any attributes/fieldnames, otherwise it will not work.

        f_out = csv.DictWriter(f_out,fieldnames=fieldnames) # use the DictWriter module and declare your fieldnames
        f_out.writeheader() # write the header
        for elem in transposons: # iterate over your structure (transposons) and for every element write a row
            f_out.writerow(elem.__dict__) # elem.__dict__ is how we access all of the attributes of an instance

def te_distribution():
    with open(TE_Table,'r') as f_in:
        reader = csv.reader(f_in, delimiter=',')
        type_dict = {}
        family_dict = {}
        next(f_in)
        for row in reader:
            the_type = str(row[5])
            the_family = str(row[6])

            if the_type in type_dict:
                type_dict[the_type] += 1
            else:
                type_dict[the_type] = 0

            if the_family in family_dict:
                family_dict[the_family] += 1
            else:
                family_dict[the_family] = 0

    type_list = []
    for key,val in type_dict.items():
        type_list.append(key)

    family_list = []
    for key,val in family_dict.items():
        family_list.append(key)

    with open(Type_Table,'w') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=type_list)
        writer.writeheader()
        writer.writerow(type_dict)

    with open(Family_Table,'w') as f_out:
        writer = csv.DictWriter(f_out,fieldnames=family_list)
        writer.writeheader()
        writer.writerow(family_dict)


class Genic_Element(object):
    def __init__(self, number, chromosome, start, stop ):
        self.number = number
        self.chromosome = chromosome
        self.start = int(start)
        self.stop = int(stop)

    def getNumber(self):
        return self.number

    def getChromosome(self):
        return self.chromosome

    def getStart(self):
        return self.start

    def getStop(self):
        return self.stop

    def getLength(self):
        return self.length

class TE(Genic_Element):
    def __init__(self, number, chromosome, start, stop, length, te_type, family):
        super().__init__(number, chromosome, start, stop)
        self.te_type = te_type
        self.family = family
        self.length = length

    def getTe_Type(self):
        return self.te_type

    def getFamily(self):
        return self.family

def run_all():
    te_handler()
    te_distribution()

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

if __name__ == '__main__':
    run_all()
