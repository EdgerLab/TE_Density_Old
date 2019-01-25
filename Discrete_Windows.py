# By Scott Teresi
# =============================================
import sys
import pandas as pd
import numpy as np

def discrete(current_density, current_window, previous_density, step = 500):
    discrete_density = (current_density * current_window) / (current_window - step) - previous_density
    return discrete_density

def set_up():
    file_name1 = sys.argv[1]
    file_name2 = sys.argv[2]
    df1=pd.read_csv(file_name1, sep=',', delimiter = None,header='infer', engine='python')
    df2=pd.read_csv(file_name2, sep=',', delimiter= None, header='infer', engine='python')

    #df1=pd.read_csv('500_xx00_density_data.csv', sep=',', delimiter = None, header='infer', engine='python')
    #df2=pd.read_csv('1000_xx00_density_data.csv', sep=',', delimiter = None, header='infer', engine='python')

    df1.replace('', np.NaN)
    df2.replace('', np.NaN)
    df1_length = len(df1.index)
    df2_length = len(df2.index)

    TEs = list(df1)
    list_to_remove = ['maker_name','chromosome','start','stop','length','prox_left','prox_right','they_are_inside','number','window_size']
    for item in list_to_remove:
        TEs.remove(item)

    if df1_length != df2_length:
        raise ValueError

    for i in range(0,df1_length):
        for item in TEs:
            new_col_name = 'Change_From_Last_Window '+ item
            current_density = np.float(df2.at[i,item])
            current_window = np.int(df2.at[i,'window_size'])
            previous_density = np.float(df1.at[i,item])

            difference = discrete(current_density, current_window, previous_density)
            df2.at[i,new_col_name] = difference

    if df1.at[0,'window_size'] == 500:
        for i in range(0,df1_length):
            for item in TEs:
                new_col_name = 'Change_From_Last_Window '+ item
                df1.at[i,new_col_name] = np.nan

    df1.to_csv('Discrete_' + file_name1)
    df2.to_csv('Discrete_' + file_name2)

set_up()
