import re, os, time, csv
files = ['500_xx00_density_data.csv','500_xx01_density_data.csv','500_xx02_density_data.csv','500_xx03_density_data.csv','500_xx04_density_data.csv']
updated_data = 'cleaned_data.csv'
dlm = ','

def cleaner():
    for item in files:
        moving_filename = item.split('_')[0] + '_' + updated_data # creates an appropriate filename for the output
        try:
             os.remove(moving_filename)
        except FileNotFoundError:
            pass

        with open(item, 'r') as f_in:
            first = True
            for row in f_in:
                row = row.rstrip()
                if first:
                    row = re.split(',',row)
                    row[len(row)-1] = row[len(row)-1].replace('}','')
                    row[0] = row[0].replace('{','')

                    my_list = []
                    for col in row:
                        col = col.split(': ')
                        col = col[0]
                        my_list.append(col)
                        first = False

                    with open(moving_filename, 'a') as f_out:
                        writer = csv.writer(f_out, delimiter=',')
                        writer.writerow(my_list)
                else:
                    row = re.split(',',row)
                    row[len(row)-1] = row[len(row)-1].replace('}','')
                    row[0] = row[0].replace('{','')
                    my_list = []

                    for col in row:
                        col = col.split(': ')
                        col = col[1]
                        my_list.append(col)

                    with open(moving_filename, 'a') as f_out:
                        writer = csv.writer(f_out, delimiter=',')
                        writer.writerow(my_list)

if __name__ == '__main__':
    cleaner()





