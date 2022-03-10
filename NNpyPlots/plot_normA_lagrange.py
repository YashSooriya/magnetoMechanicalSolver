import plots
import numpy as np
from normA import NormA


def main():
    Options = {'POD': 'PODI','neural_network': True, 'feed_forward_net': True, 'segmented': False, 'per_mode': True, 'N_s': 180, 'm': 20,
            'N_o': 40, 'x_axis_logged': False, 'No_segments': 20, 'layers': 1, 'chosen_sample': 'log space',
            'neurons': 1, 'solver': 'adam', 'activation': 'identity', 'ffn_layers': 2, 'ffn_neurons': 16, 'ffn_solver': 'trainbr',
               'custom_saveload':True, 'save_name': '', 'load_name': '', 'min_grad': 1e-9}
    
    # freq_out = np.linspace(15,15+(Options['N_o']-1)*10, Options['N_o'])
    # freq_out = np.linspace(5, 105, N_o)
    freq_out = np.linspace(15, 5000, Options['N_o'])
    # freq_out = np.logspace(np.log10(5),np.log10(5000),Options['N_s'])
    # freq_out = np.log10(np.linspace(15, 5000, Options['N_o']))
    Options['freqout'] = freq_out


    globals().update(Options)
    min_grad_n = -int(np.log10(min_grad))

    Options['load_name'] = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}_freqlog_mingrad{min_grad_n}.mat'


    Options['save_name'] = f'normA/NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog_mingrad{min_grad_n}'
    normA = NormA(Options)
    y_NN, save_name_NN = normA.load()
    plot_options = normA.plot_options()
    labels_nn = plot_options['labels']
    xlabel = plot_options['xlabel']
    ylabel = plot_options['ylabel']

    Options['neural_network'] = False
    Options['load_name'] = f'FrequencySweepMHIGradXNormA_lagrange_15to5000Hzlogspace_damp2e-3_q3_p3_Ns{N_s}PODI_serial_OVCNS_No40_m{m}.mat'
    normA_lagr = NormA(Options)
    y_lagr, save_name_lagr = normA_lagr.load()
    plot_options_lagr = normA_lagr.plot_options()
    labels_lagr = plot_options_lagr['labels']

    y = y_lagr + y_NN
    labels = labels_lagr + labels_nn

    save_name = Options['save_name']
    plots.list_log_plot(freq_out, y, xlabel, ylabel, labels, save_name, lgndlocation=1)




if __name__ == '__main__':
    main()
