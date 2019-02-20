# By Scott Teresi
# =======================================================
import sys
import pandas as pd
import numpy as np
import os
import multiprocessing
from multiprocessing import Process
from threading import Thread
import time
from collections import defaultdict

def discrete(current_density, current_window, previous_density, step = 500):
    '''Mathematical function to compute the change in density between windows'''
    discrete_density = (current_density * current_window) / (current_window - step) - previous_density # formula for the density calculation to remove previous window
    return discrete_density # return value

def set_up(file_name1, file_name2):
    df1=pd.read_csv(file_name1, sep=',', delimiter = None, header='infer', engine='python') # read the data
    df2=pd.read_csv(file_name2, sep=',', delimiter = None, header='infer', engine='python')


    df1.replace('', np.NaN) # for objects that have no value, replace with NULL value in the PandaFrame
    df2.replace('', np.NaN)
    df1_length = len(df1.index) # get length for later iteration
    df2_length = len(df2.index)

    if df1_length != df2_length: # value check to makemsure PandaFrames are the same length. This is expected
        raise ValueError

    TEs = list(df1) # generate fieldnames list from PandaFrame
    list_to_remove = ['maker_name','chromosome','start','stop','length','prox_left','prox_right','they_are_inside','number','window_size']
        # fieldnames to remove from fieldnames list

    for item in list_to_remove:
        TEs.remove(item) # we will only run the calculations on the filtered fieldnames


    for i in range(0,df1_length): # iterate over range of numbers to correspond to PandaFrame
        for item in TEs: # iterate over acceptable fieldnames
            new_col_name = 'Change_From_Last_Window '+ item # generate fieldname for new data column
            current_density = np.float(df2.at[i,item]) # declare as Numpy object, density equal to fieldname (item) density at ith row


            current_window = np.int(df2.at[i,'window_size'])  # declare as Numpy object, window equal to fieldname (item) window size at ith row (same for whole file)

                #  TO DO:
                    #Could streamline to declare once


            previous_density = np.float(df1.at[i,item]) # declare as Numpy object, previous density equal to item density at ith row

            difference = discrete(current_density, current_window, previous_density) # calculate the difference from previous window and assign it to difference
            df2.at[i,new_col_name] = difference # Add to df2 a new column, at the ith row, the value of difference

    if df1.at[0,'window_size'] == 500: # Add the columns to 500 bp files, just to have the columns, and fill the values with NULL
        for i in range(0,df1_length):
            for item in TEs:
                new_col_name = 'Change_From_Last_Window '+ item
                df1.at[i,new_col_name] = np.nan





    # TO DO
        # drop the unneccessary columns since I only want the discrete jumps, to save space and RAM when computing
    #list_to_remove_2 = list_to_remove.remove('maker_name')
#['maker_name','chromosome','start','stop','length','prox_left','prox_right','they_are_inside','number','window_size']

    #df1.drop(columns = ['chromosome','start','stop','length','prox_left','prox_right','they_are_inside','number','window_size'])
    #df2.drop(columns = ['chromosome','start','stop','length','prox_left','prox_right','they_are_inside','number','window_size'])

    df1.to_csv('Discrete_' + file_name1) # write the PandaFrames to a csv
    df2.to_csv('Discrete_' + file_name2)

def show_status(status_queue):
    """Print to screen the percent completion of each response."""

    time.sleep(10)  # MAGIC NUMBER wait for printouts from the process start
    status_map = defaultdict(tuple)
    while True:  # MAGIC NUMBER make sure to use a daemon if no kill event
        response = status_queue.get(block=True)
        try:
            percent, my_id = response
        except ValueError as ve:
            msg = ("show_status expecting 3 values in response but got {}"
                   .format(response))
            logger.debug(msg)
            continue

        status_map[str(my_id)] = percent
        status = "percent complete: "
        for id, complete in status_map.items():
            status += id + ": {:.2f}".format(complete) + ",  "

        sys.stdout.write('\r')
        sys.stdout.write(status)
        sys.stdout.flush()

def run_all():
    os.chdir('CAMDATA')
    my_inputs = ['xx00','xx01','xx02','xx03','xx04','xx05','xx06'] # determine the chromosome headings of files to loop over
    count = 500 # the jump in windows of the file scheme
    my_dict = {}

    for item in my_inputs:
        last = str(count+500) + '_' + item + '_density_data.csv'
        while count <= 10000:

            my_file = str(count) + '_' + item + '_density_data.csv'
            my_dict[last] = str(count) + '_' + item + '_density_data.csv'
            count += 500
            last = my_file
        count = 500


    file_list = []
    for item in my_inputs:
        file_list.append(item[0])

    status_queue = multiprocessing.Queue()
    status_thread = Thread(target=show_status, args=(status_queue,))
    status_thread.daemon = True

    p_list = []
    for f_name in my_dict.keys():
        p = Process(target=set_up, args=(f_name, my_dict[f_name],))
        p.start()
        p_list.append(p)

    status_thread.start()

    for p in p_list:
        current_pid = p.pid
        p.join(None)
        print("pid '{}' joined!".format(current_pid))

run_all()
