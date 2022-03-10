import numpy as np
from scipy.io import loadmat
import csv
from pathlib import Path
from pandas import DataFrame, read_csv
from normA import NormA

def main():
    Options = {'neural_network': False, 'nftool': True, 'segmented': True, 'per_mode': False, 'N_s': 599, 'm': 20,
           'N_o': 40, 'x_axis_logged': False, 'No_segments': 20, 'layers': 2,
               'neurons': 2, 'solver': 'adam', 'activation': 'identity', 'nftool_neurons': 50, 'nftool_solver': 'trainbr'}
    # freq_out = np.linspace(15,15+(N_o-1)*10, N_o)
    # freq_out = np.linspace(5, 105, N_o)
    freq_out = np.linspace(15, 5000, Options['N_o'])
    Options['freqout'] = freq_out
    
    normA_lagr = NormA(Options)
    y_lagr, save_name_lagr = normA_lagr.load()
    
    Options['neural_network'] = True
    normA_nn = NormA(Options)
    y_nn, save_name_nn = normA_nn.load()
    
    table_entry = calc_diff(y_lagr, y_nn)
    table_entry.insert(0, save_name_nn)

    file_name = 'normA_differences.csv'
    write_to_file(table_entry, file_name)
    print("\n")
    display_csv(file_name)
    
    print("---------------------------------------------------------------")
    print(table_entry)
    print("---------------------------------------------------------------")
    print(f"Successfully calculated difference in norm A for: \n {Options}")


def display_csv(file_name):
    df = read_csv(file_name, sep=',')
    df.drop_duplicates(subset=None, inplace=True)
    print(df)
    df.to_csv(file_name, index=False)
    
    
def write_to_file(entry_list, file_name):
    fields = ['NN_options', '4K', '77K', 'OVC']
    if not Path(file_name).is_file():
        with open(file_name, 'w') as filehandle:
            write = csv.writer(filehandle)
            write.writerow(fields)
    
    with open(file_name, 'a') as filehandle:
        write = csv.writer(filehandle)
        write.writerow(entry_list)
    
def calc_diff(A1, A2):
    table_entry = []
    for i in range(3):
        diff = np.nansum((A1[i] - A2[i]) ** 2) ** 0.5
        table_entry.append(diff / (np.nansum(A1[i] ** 2) ** 0.5))
    return table_entry

if __name__ == '__main__':
    main()
