import funcs
import plots
import numpy as np
from scipy.io import loadmat
from normA import NormA

def main():
    Options = {'POD': 'PODI','neural_network': True, 'feed_forward_net': True, 'segmented': False, 'per_mode': True, 'N_s': 23, 'm': 20,
            'N_o': 40, 'x_axis_logged': False, 'No_segments': 20, 'layers': 1, 'chosen_sample': 'log space',
            'neurons': 1, 'solver': 'adam', 'activation': 'identity', 'ffn_layers': 2, 'ffn_neurons': 32, 'ffn_solver': 'trainbr',
               'custom_saveload':True,'save_name': '', 'load_name': '', 'min_grad': 1e-10}
    
    # freq_out = np.linspace(15,15+(Options['N_o']-1)*10, Options['N_o'])
    # freq_out = np.linspace(5, 105, N_o)
    freq_out = np.linspace(15, 5000, Options['N_o'])
    # freq_out = np.logspace(np.log10(5),np.log10(5000),Options['N_s'])
    # freq_out = np.log10(np.linspace(15, 5000, Options['N_o']))
    Options['freqout'] = freq_out


    REG = False
    MIN_GRAD = False
    LAGRANGE = True



    globals().update(Options)

    min_grad_n = -int(np.log10(min_grad))


    Options['load_name'] = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}_freqlog_mingrad{min_grad_n}.mat'
    
    
    normA = NormA(Options)
    y1, save_name_1 = normA.load()
    plot_options = normA.plot_options()
    xlabel = plot_options['xlabel']
    ylabel = plot_options['ylabel']
    
    Options['load_name'] = f'FrequencySweepMHIGradXNormA_approxD_5to5000Hzlogspace_damp2e-3_q3_p3_Ns{N_s}PODI_serial_OVCNS_NsNo_m20.mat'
    normA_approxD = NormA(Options)
    y_approxD, _ = normA_approxD.load()
    
    approxD_freqs = np.logspace(np.log10(5),np.log10(5000),N_s)

    if REG:
        save_name = f'normA/NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog_approxDcompare_regs'
        
        Options['load_name'] = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}_freqlog_reg10.mat'
        normA = NormA(Options)
        y2, save_name_2 = normA.load()
    
        Options['load_name'] = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}_freqlog_reg20.mat'
        normA = NormA(Options)
        y3, save_name_3 = normA.load()
    
        labels = ['approx D', 'NN 0% Regularization', 'NN 10% Regularization', 'NN 20% Regularization']
        ys = [y_approxD, y1, y2, y3]
        
    elif MIN_GRAD:
        save_name = f'normA/NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog_approxDcompare_mingrads'
        Options['load_name'] = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}_freqlog_mingrad8.mat'
        normA = NormA(Options)
        y2, save_name_2 = normA.load()

        Options['load_name'] = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}_freqlog_mingrad9.mat'
        normA = NormA(Options)
        y3, save_name_3 = normA.load()
        
        labels = ['approx D', 'NN 10e-7 minimum grad', 'NN 10e-8 minimum grad', 'NN 10e-9 minimum grad']
        ys = [y_approxD, y1, y2, y3]
        
    elif LAGRANGE:
        save_name = f'normA/NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog_approxDcompare_lagrange'
        Options['neural_network'] = False
        Options['load_name'] = f'FrequencySweepMHIGradXNormA_lagrange_15to5000Hzlogspace_damp2e-3_q3_p3_Ns{N_s}PODI_serial_OVCNS_No40_m{m}.mat'
        normA = NormA(Options)
        y_lagr, save_name_2 = normA.load()
        labels = ['approx D', 'Neural Network', 'Lagrange']
        ys = [y_approxD, y1, y_lagr]


    shield_4k = [i[0] for i in ys]
    
    plots.list_scatter_line_log_multiple_plot(freq_out,approxD_freqs,shield_4k, xlabel, ylabel, labels, save_name)
    

if __name__ == '__main__':
    main()
