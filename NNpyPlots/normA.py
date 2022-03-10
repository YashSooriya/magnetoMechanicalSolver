import numpy as np
from scipy.io import loadmat

class NormA:

    def __init__(self, Options):
        self.Options = Options

    def load(self):
        if self.Options['custom_saveload']==True:
            load_name = self.Options['load_name']
            save_name = self.Options['save_name']
            load_directory = 'data/normA/'

        else:
            load_name, save_name = self.__file_names()
            load_directory = '../data/normA/old/'


        load = load_directory + load_name
        save_directory = 'normA/'
        save = save_directory + save_name

        data = loadmat(load)
        print(f'Loaded {load}')
        struct = data['IntegratedNormA']
        output = [struct[0, 0][i] for i in range(3)]

        return output, save

    def plot_options(self):
        
        if self.Options['neural_network']:
            labels_list = ['4K Neural Network','77K Neural Network', 'OVC Neural Network']
        else:
            if self.Options['POD'] == 'PODI':
                labels_list = ['4K Lagrange','77K Lagrange', 'OVC Lagrange']
            else:
                labels_list = ['4K PODP','77K PODP', 'OVC PODP']
                
            
        x_label = 'Frequency (Hz)'
        y_label = '$||\mathbf{A}^{AC}_{\epsilon, hp}||$ (Vs$m^-1$)'
        options_dict = {'xlabel': x_label, 'ylabel': y_label, 'labels':labels_list}

        return options_dict

    def __file_names(self):
        neural_network = self.Options['neural_network']
        feed_forward_net = self.Options['feed_forward_net']
        ffn_layers = self.Options['ffn_layers']
        ffn_neurons = self.Options['ffn_neurons']
        ffn_solver = self.Options['ffn_solver']
        segmented = self.Options['segmented']
        per_mode = self.Options['per_mode']
        N_s = self.Options['N_s']
        m = self.Options['m']
        N_o = self.Options['N_o']
        x_axis_logged = self.Options['x_axis_logged']
        freqout = self.Options['freqout']
        No_segments = self.Options['No_segments']
        layers = self.Options['layers']
        neurons = self.Options['neurons']
        solver = self.Options['solver']
        activation = self.Options['activation']
        POD = self.Options['POD']
        chosen_sample = self.Options['chosen_sample']
        
        if neural_network:
            if segmented:
                custom_filename = f'FrequencySweepMHIGradXNormA_NN{No_segments}Segments_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}.mat'
            elif per_mode:
                custom_filename = f'FrequencySweepMHIGradXNormA_NN{m}Modes_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}.mat'
                if x_axis_logged:
                    custom_filename = f'FrequencySweepMHIGradXNormA_NN{m}Modes_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_logged_OVCNS_No{N_o}.mat'
            else:
                custom_filename = f'FrequencySweepMHIGradXNormA_NN_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}.mat'
            if feed_forward_net:
                if per_mode:
                    if chosen_sample == "log space":
                        custom_filename = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hzlogspace_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}.mat'
                    elif chosen_sample == 'marcos':
                        custom_filename = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}.mat'                        
                else:
                    custom_filename = f'FrequencySweepMHIGradXNormA_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_allModes_No{N_o}.mat'
        else:     
            custom_filename = f'FrequencySweepMHIGradXNormA_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}_m{m}.mat'

        if neural_network:
            if segmented:
                custom_savename = f'NormA_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_{No_segments}segmented_{int(freqout[0])}to{int(freqout[-1])}Hz'
            if per_mode:
                custom_savename = f'NormA_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_per{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz'
                if x_axis_logged:
                    custom_savename = f'NormA_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_per{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_xlogged'
            if feed_forward_net:
                if per_mode:
                    custom_savename = f'NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz'
                else:
                    custom_savename = f'NormA_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_allModes'

            else:
                custom_savename = f'NormA_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_{int(freqout[0])}to{int(freqout[-1])}Hz'

        else:
            custom_savename = f'NormA_{POD}_Ns{N_s}_No{N_o}_{int(freqout[0])}to{int(freqout[-1])}Hz_m{m}'
            if x_axis_logged:
                custom_savename = f'NormA_{POD}_Ns{N_s}_No{N_o}_{int(freqout[0])}to{int(freqout[-1])}Hz_xlogged'

        return custom_filename, custom_savename

        
