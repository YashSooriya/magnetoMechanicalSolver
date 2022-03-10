import numpy as np
from scipy.io import loadmat

class Energy:

    def __init__(self, Options):
        self.Options = Options

    def load(self):
        load_name, save_name = self.__file_names()

        load_directory = '../data/powerEnergy/'
        load = load_directory + load_name
        save_directory = 'energy/'
        save = save_directory + save_name

        # density
        rho4K = 2710
        rho77K = 2698
        rhoOVC = 7900

        omega = 31415.93


        data = loadmat(load)
        print(f'Loaded {load}')
        struct = data['IntegratedFields'][0, 0]
        
        out4K = np.square(struct['DisplacementNorm4K']) * 2 * rho4K * (omega ** 2)
        out77K = np.square(struct['DisplacementNorm77K']) * 2 * rho77K * (omega ** 2)
        outOVC = np.square(struct['DisplacementNormOVC']) * 2 * rhoOVC * (omega ** 2)
        output = [out4K, out77K, outOVC]
        
        return output, save

    def plot_options(self):
        
        if self.Options['neural_network']:
            labels_list = ['4K Neural Network','77K Neural Network', 'OVC Neural Network']
        else:
            labels_list = ['4K Lagrange','77K Lagrange', 'OVC Lagrange']
            
        x_label = 'Frequency (Hz)'
        y_label = 'Kinetic energy (J)'
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
        
        if neural_network:
            if segmented:
                custom_filename = f'FrequencySweepMHIGradXPowerEnergy_NN{No_segments}Segments_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}.mat'
            elif per_mode:
                custom_filename = f'FrequencySweepMHIGradXPowerEnergy_NN{m}Modes_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}.mat'
                if x_axis_logged:
                    custom_filename = f'FrequencySweepMHIGradXPowerEnergy_NN{m}Modes_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_logged_OVCNS_No{N_o}.mat'
            else:
                custom_filename = f'FrequencySweepMHIGradXPowerEnergy_NN_{solver}_{activation}_l{layers}_n{neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}.mat'
            if feed_forward_net:
                if per_mode:
                    custom_filename = f'FrequencySweepMHIGradXPowerEnergy_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}.mat'
                else:
                    custom_filename = f'FrequencySweepMHIGradXPowerEnergy_NN_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_allModes_No{N_o}.mat'
        else:     
            custom_filename = f'FrequencySweepMHIGradXPowerEnergy_{int(freqout[0])}to{int(freqout[-1])}Hz_damp2e-3_q3_p3_Ns{N_s}{POD}_serial_OVCNS_No{N_o}.mat'

        if neural_network:
            if segmented:
                custom_savename = f'Energy_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_{No_segments}segmented_{int(freqout[0])}to{int(freqout[-1])}Hz'
            if per_mode:
                custom_savename = f'Energy_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_per{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz'
                if x_axis_logged:
                    custom_savename = f'Energy_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_per{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_xlogged'
            if feed_forward_net:
                if per_mode:
                    custom_savename = f'Energy_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz'
                else:
                    custom_savename = f'Energy_{POD}_ffn_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_No{N_o}_{m}modes_{int(freqout[0])}to{int(freqout[-1])}Hz_allModes'

            else:
                custom_savename = f'Energy_{POD}_{solver}_{activation}_l{layers}_n{neurons}_Ns{N_s}_No{N_o}_{int(freqout[0])}to{int(freqout[-1])}Hz'

        else:
            custom_savename = f'Energy_{POD}_Ns{N_s}_No{N_o}_{int(freqout[0])}to{int(freqout[-1])}Hz'
            if x_axis_logged:
                custom_savename = f'Energy_{POD}_Ns{N_s}_No{N_o}_{int(freqout[0])}to{int(freqout[-1])}Hz_xlogged'

        return custom_filename, custom_savename

        
