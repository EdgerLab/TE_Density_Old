# By Scott Teresi --- July 2018
#---------------------------------------------
import csv
import glob

# This function groups the TE's of my file into bins for GO enrichment.
# I use this file to partition the density of my genes

file_list = glob.glob("/home/scott/Documents/Uni/Research/transposon_2.0/Camarosa_TE/CAMDATA/*density_data.csv")

#-----------------------
# if fieldnames is not supplied to DictReader then it uses the first line to make the fieldnames

def test():
    #for item in file_list:


    with open('GO_Groupings.csv') as f_out:
        writer = csv.DictWriter(f_out,fieldnames = fieldnames)
        with open("TE_Data.csv") as csv_file:
            reader = csv.DictReader(csv_file)
            for row in reader: # refer to the variables through string name
                name = row['maker_name']
                window_size = int(row['window_size'])
                classification = row['variable']
                value = float(row['value'])
                combined_name = str(window_size)+ '_' + classification

                if value < Lower[combined_name][2] or value > Upper[combined_name][2]:
                    # Keep the data






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




if __name__ == '__main__':
    whitelist_maker()
    test()
