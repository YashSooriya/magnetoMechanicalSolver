import funcs
import plots
import numpy as np
from scipy.io import loadmat
from energy import Energy

def main():
    Options = {'POD': 'PODP','neural_network': True, 'feed_forward_net': True, 'segmented': False, 'per_mode': True, 'N_s': 180OB, 'm': 20,
           'N_o': 499, 'x_axis_logged': False, 'No_segments': 20, 'layers': 1,
               'neurons': 1, 'solver': 'adam', 'activation': 'identity', 'ffn_layers': 1, 'ffn_neurons': 8, 'ffn_solver': 'trainbr'}

    LAGRANGE_PODI = False
    NEURAL_NETWORK = True
    
    freq_out = np.linspace(15,15+(Options['N_o']-1)*10, Options['N_o'])
    # freq_out = np.linspace(5, 105, N_o)
    # freq_out = np.linspace(15, 5000, Options['N_o'])
    Options['freqout'] = freq_out

    energy_nn = Energy(Options)
    y_nn, save_name_nn = energy_nn.load()
    plot_options_nn = energy_nn.plot_options()
    labels_nn = plot_options_nn['labels']
    

    xlabel = plot_options_nn['xlabel']
    ylabel = plot_options_nn['ylabel']


    if LAGRANGE_PODI:
        Options['neural_network'] = False
        
        energy_lagr = Energy(Options)
        y_lagr, save_name_lagr = energy_lagr.load()
        plot_options_lagr = energy_lagr.plot_options()
        labels_lagr = plot_options_lagr['labels']
        
    
        labels = labels_lagr + labels_nn
        y = y_lagr + y_nn
        save_name = save_name_nn + '_lagrCompare'


    
    if LAGRANGE_PODI and NEURAL_NETWORK:
        plots.list_log_plot(freq_out, y, xlabel, ylabel, labels, save_name, lgndlocation=0)
#        plots.list_log_plot(freq_out, difference, xlabel, ylabel, diff_labels, save_name_diff)
    elif LAGRANGE_PODI:
        plots.list_log_plot(freq_out, y_lagr, xlabel, ylabel, labels_lagr, save_name_lagr)
    elif NEURAL_NETWORK:
        plots.list_log_plot(freq_out, y_nn, xlabel, ylabel, labels_nn, save_name_nn)


    print(f"Successfully plotted for: \n {Options}")


if __name__ == '__main__':
    main()
