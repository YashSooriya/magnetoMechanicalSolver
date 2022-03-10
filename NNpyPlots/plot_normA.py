import plots
import numpy as np
from normA import NormA


def main():    
    Options = {'POD': 'PODI', 'neural_network': True, 'feed_forward_net': True, 'segmented': False, 'per_mode': True, 'N_s': 90, 'm': 20,
               'N_o': 40, 'x_axis_logged': False, 'No_segments': 20, 'layers': 1, 'chosen_sample': 'log space',
               'neurons': 1, 'solver': 'adam', 'activation': 'identity', 'ffn_layers': 2, 'ffn_neurons': 16, 'ffn_solver': 'trainbr',
               'custom_saveload': True, 'save_name': '', 'load_name': ''}

    NO_COMPARISON = True
    LAGRANGE_OR_PODP = False
    NEURAL_NETWORK = False
    DIFFERENCE = False
    LAGRANGE_OR_PODP_M_PER_SHIELD = False
    APPROX_D_COMPARE = False
    extra_save = ''
    

    # freq_out = np.linspace(15,15+(Options['N_o']-1)*10, Options['N_o'])
    # freq_out = np.linspace(5, 105, N_o)
    freq_out = np.linspace(15, 5000, Options['N_o'])
    # freq_out = np.logspace(np.log10(5),np.log10(5000),Options['N_s'])
    # freq_out = np.log10(np.linspace(15, 5000, Options['N_o']))
    Options['freqout'] = freq_out

    neural_network = Options['neural_network']
    feed_forward_net = Options['feed_forward_net']
    ffn_layers = Options['ffn_layers']
    ffn_neurons = Options['ffn_neurons']
    ffn_solver = Options['ffn_solver']
    segmented = Options['segmented']
    per_mode = Options['per_mode']
    N_s = Options['N_s']
    m = Options['m']
    N_o = Options['N_o']
    x_axis_logged = Options['x_axis_logged']
    freqout = Options['freqout']
    No_segments = Options['No_segments']
    layers = Options['layers']
    neurons = Options['neurons']
    solver = Options['solver']
    activation = Options['activation']
    POD = Options['POD']
    chosen_sample = Options['chosen_sample']

    if neural_network:
        Options['load_name'] = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}_freqlog.mat'
    else:
        Options['load_name'] = f'FrequencySweepMHIGradXNormA_lagrange_15to5000Hzlogspace_damp2e-3_q3_p3_Ns{N_s}PODI_serial_OVCNS_No{N_o}_m{m}.mat'
    # Options['load_name'] = f'FrequencySweepMHIGradXNormA_approxD_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_NsNo_m{m}.mat'
    # Options['load_name'] = f'FrequencySweepMHIGradXNormA_approxD_5to5000Hzlogspace_damp2e-3_q3_p3_Ns180PODI_serial_OVCNS_NsNo_m20.mat'

    Options['save_name'] = f'NormA_{POD}_approxD_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz'

    # Options['save_name'] = f'NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog'

    if APPROX_D_COMPARE:
        normA = NormA(Options)
        y, save_name_1 = normA.load()
        plot_options = normA.plot_options()
        xlabel = plot_options['xlabel']
        ylabel = plot_options['ylabel']
        
        if Options['custom_saveload'] == True:
            save_name = f'normA/NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog_approxDcompare'
        else:
            save_name = save_name_1 + '_approxDcompare'

        Options['custom_saveload'] = True
        Options['load_name'] = f'FrequencySweepMHIGradXNormA_approxD_5to5000Hzlogspace_damp2e-3_q3_p3_Ns{N_s}PODI_serial_OVCNS_NsNo_m20.mat'
        normA_approxD = NormA(Options)
        y_approxD, _ = normA_approxD.load()

        x2 = np.logspace(np.log10(5),np.log10(5000),Options['N_s'])

        ys = list(zip(y, y_approxD))

        if Options['neural_network'] == True:
            labels = [('4K NN', '4K approx D'), ('77K NN', '77K approx D'), ('OVC NN', 'OVC approx D')]
            save_name = f'normA/NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog_approxDcompare'

        else:
            labels = [('4K Lagrange', '4K approx D'), ('77K Lagrange', '77K approx D'), ('OVC Lagrange', 'OVC approx D')]
            save_name = f'normA/NormA_{POD}_lagrange_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_freqlog_approxDcompare'



        plots.list_scatter_line_log_plot(freq_out, x2, ys, xlabel, ylabel, labels, save_name, marker_size=4)

    if LAGRANGE_OR_PODP_M_PER_SHIELD:
        m_values = [5, 10, 15, 20]
        y_ms = []
        Options['neural_network'] = False
        xlabel = ''
        ylabel = ''
        for i in range(len(m_values)):
            Options['m'] = m_values[i]
            normA_m = NormA(Options)
            y_m, _ = normA_m.load()
            plot_options_m = normA_m.plot_options()

            y_ms.append(y_m)
            xlabel = plot_options_m['xlabel']
            ylabel = plot_options_m['ylabel']
        

        y_m_4k, y_m_77k, y_m_ovc = lagr_m_per_shield(y_ms)
        labels = [f'm = {i}' for i in m_values]
        save_name = f'normA/NormA_Ns{Options["N_s"]}_No{Options["N_o"]}_{int(freq_out[0])}to{int(freq_out[-1])}Hz_varyModes'
        plots.list_log_plot(freq_out, y_m_4k ,xlabel, ylabel, labels, save_name + '_4k')
        plots.list_log_plot(freq_out, y_m_77k ,xlabel, ylabel, labels, save_name + '_77k')
        plots.list_log_plot(freq_out, y_m_ovc ,xlabel, ylabel, labels, save_name + '_ovc')

        

                
    else:
        if NO_COMPARISON:
            normA = NormA(Options)
            y, save_name = normA.load()
            plot_options = normA.plot_options()
            labels = plot_options['labels']
            xlabel = plot_options['xlabel']
            ylabel = plot_options['ylabel']
            plots.list_log_plot(freq_out, y, xlabel, ylabel, labels, save_name, lgndlocation=1)

        if NEURAL_NETWORK:
            Options['neural_network'] = True
            normA_nn = NormA(Options)
            y_nn, save_name_nn = normA_nn.load()
            plot_options_nn = normA_nn.plot_options()
            labels_nn = plot_options_nn['labels']
            xlabel = plot_options_nn['xlabel']
            ylabel = plot_options_nn['ylabel']

        
        if LAGRANGE_OR_PODP:
            Options['neural_network'] = False
            Options['custom_saveload'] = False
            Options['m'] = 20
            normA_lagr = NormA(Options)
            y_lagr, save_name_lagr = normA_lagr.load()
            plot_options_lagr = normA_lagr.plot_options()
            labels_lagr = plot_options_lagr['labels']
        
        
            xlabel = plot_options_lagr['xlabel']
            ylabel = plot_options_lagr['ylabel']
    
            labels = labels_lagr + labels_nn
            y = y_lagr + y_nn
            save_name = save_name_nn + extra_save
    
        if LAGRANGE_OR_PODP and NEURAL_NETWORK:
            if DIFFERENCE:
                difference = lagr_nn_difference(y_lagr, y_nn, average=False)
                print(difference)
                save_name_diff = save_name + '_diff'
                diff_labels = ['4K Difference', '77K Difference', 'OVC Difference']
                plots.list_log_plot(freq_out, difference,xlabel, ylabel, diff_labels, save_name_diff)
            else:
                save_name = save_name + '_lagrCompare'
                plots.list_log_plot(freq_out, y, xlabel, ylabel, labels, save_name, lgndlocation=1)
                #        plots.list_log_plot(freq_out, difference, xlabel, ylabel, diff_labels, save_name_diff)
        elif LAGRANGE_OR_PODP:
            plots.list_log_plot(freq_out, y_lagr, xlabel, ylabel, labels_lagr, save_name_lagr)
        elif NEURAL_NETWORK:
            plots.list_log_plot(freq_out, y_nn, xlabel, ylabel, labels_nn, save_name_nn)
        
        print(f"Successfully plotted for: \n {Options}")


def lagr_nn_difference(y1, y2, average=True):
    norm_diff_shields = np.zeros((3, y1[0].shape[0]))
    for i in range(3):
        diff = (y1[i] - y2[i]) / y1[i]
        norm_diff_shields[i, :] = diff.ravel() 
            

    if average:
        ret = np.nanmean(norm_diff_shields, axis=0)
    else:
        ret = norm_diff_shields
    return ret

    
def calc_diff(A1, A2):
    table_entry = []
    for i in range(3):
        diff = np.nansum((A1[i] - A2[i]) ** 2) ** 0.5
        table_entry.append(diff / (np.nansum(A1[i] ** 2) ** 0.5))
    return table_entry


def lagr_m_per_shield(y_ms):
    y_4k = []
    y_77k = []
    y_ovc = []
    for j in range(len(y_ms)):
        y_4k.append(y_ms[j][0])
        y_77k.append(y_ms[j][0])
        y_ovc.append(y_ms[j][0])
                
    return y_4k, y_77k, y_ovc
    
    

if __name__ == '__main__':
    main()


