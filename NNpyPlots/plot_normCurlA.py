import funcs
import plots
import numpy as np
from scipy.io import loadmat

from normCurlA import NormCurlA
def main():
    Options = {'neural_network': True, 'feed_forward_net': True, 'segmented': True, 'per_mode': True, 'N_s': 359, 'm': 20,
           'N_o': 40, 'x_axis_logged': False, 'No_segments': 20, 'layers': 2,
               'neurons': 2, 'solver': 'adam', 'activation': 'identity', 'ffn_layers': 3, 'ffn_neurons': 16, 'ffn_solver': 'trainbr'}

    LAGRANGE = True
    NEURAL_NETWORK = True
    
    # freq_out = np.linspace(15,15+(N_o-1)*10, N_o)
    # freq_out = np.linspace(5, 105, N_o)
    freq_out = np.linspace(15, 5000, Options['N_o'])
    Options['freqout'] = freq_out

    normA_nn = NormCurlA(Options)
    y_nn, save_name_nn = normA_nn.load()
    plot_options_nn = normA_nn.plot_options()

    labels_nn = plot_options_nn['labels']
    xlabel = plot_options_nn['xlabel']
    ylabel = plot_options_nn['ylabel']

    Options['neural_network'] = False
    
    normA_lagr = NormCurlA(Options)
    y_lagr, save_name_lagr = normA_lagr.load()
    plot_options_lagr = normA_lagr.plot_options()
    labels_lagr = plot_options_lagr['labels']

    
    labels = labels_lagr + labels_nn
    y = y_lagr + y_nn
    save_name = save_name_nn + '_lagrCompare'


    
    if LAGRANGE and NEURAL_NETWORK:
        plots.list_log_plot(freq_out, y, xlabel, ylabel, labels, save_name, lgndlocation=1)
#        plots.list_log_plot(freq_out, difference, xlabel, ylabel, diff_labels, save_name_diff)
    elif LAGRANGE:
        plots.list_log_plot(freq_out, y_lagr, xlabel, ylabel, labels_lagr, save_name_lagr)
    elif NEURAL_NETWORK:
        plots.list_log_plot(freq_out, y_nn, xlabel, ylabel, labels_nn, save_name_nn)
        
    
    print(f"Successfully plotted for: \n {Options}")

if __name__ == '__main__':
    main()
