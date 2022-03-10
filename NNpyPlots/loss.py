import numpy as np
from scipy.io import loadmat


class Loss:

    def __init__(self, Options):
        self.Options = Options


    def load(self):
        load_name, save_name = self.__file_names()
        
        data = loadmat(load_name)
        print(f'Loaded {load_name}')
        struct = data['tr'][0, 0]
        
        metrics = {'train_losses': struct['perf'], 'test_losses': struct['tperf'], 'epochs': struct['epoch'],
                   'stop': struct['stop'], 'best_train_loss': struct['best_perf'], 'best_test_loss': struct['best_tperf'],
                   'best_epoch': struct['best_epoch']}
        
        metrics_T = {k: np.transpose(v) for k,v in metrics.items()}

        return metrics_T, save_name

    def plot_options(self):
        x_label = 'Epochs'
        y_label = 'Loss (MSE)'
        options_dict = {'xlabel': x_label, 'ylabel': y_label}
        
        return options_dict

    def __file_names(self):
        layers = self.Options['layers']
        neurons = self.Options['neurons']
        solver = self.Options['solver']
        per_mode = self.Options['per_mode']
        N_s = self.Options['N_s']
        m = self.Options['m']
        network_no = self.Options['network_no']
        PODI = self.Options['PODI']
        PODP = self.Options['PODP']

        if per_mode:
            load_file = f'LossCurve_{solver}_l{layers}_n{neurons}_m{m}_Ns{N_s}_network{network_no}.mat'
            save_file = f'LossCurve_{solver}_l{layers}_n{neurons}_m{m}_Ns{N_s}_network{network_no}'
        else:
            load_file = f'LossCurve_{solver}_l{layers}_n{neurons}_m{m}_Ns{N_s}_allModes.mat'
            save_file = f'LossCurve_{solver}_l{layers}_n{neurons}_m{m}_Ns{N_s}_allModes'

        if PODI:
            load_directory = f'../NNdata/PODI/loss_curves/loss_curves_{solver}_l{layers}_n{neurons}_m{m}_Ns{N_s}/'
            save_directory = 'loss_curves/PODI/'

        elif PODP:
            load_directory = f'../NNdata/PODP/loss_curves/loss_curves_{solver}_l{layers}_n{neurons}_m{m}_Ns{N_s}/'
            save_directory = 'loss_curves/PODP/'

            

        load = load_directory + load_file
        save = save_directory + save_file
        return load, save
        

