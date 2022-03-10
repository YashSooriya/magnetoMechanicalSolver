import plots
import numpy as np
from scipy.io import loadmat


def main():
    Options = {'N_s': 23, 'm': 20, 'N_o': 40, 'chosen_sample': 'log space',
               'neurons': 1, 'solver': 'adam', 'activation': 'identity',
               'ffn_layers': 3, 'ffn_neurons': 16, 'ffn_solver': 'trainbr',
               'min_grad': 1e-10, 'reg': 0}

    globals().update(Options)

    freq_out = np.linspace(5, 5000, N_o)

    data_lgr = loadmat(f'data/GCurve/GCurve_lagrange_Ns{N_s}_m{m}.mat')
    G_lgr = data_lgr['VintWhole']

    data_NN = loadmat(f'data/GCurve/GCurve_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_m{m}_Ns{N_s}_reg{reg}.mat')
    G_NN = data_NN['new_G']

    G_lgr = np.absolute(G_lgr)
    G_NN = np.absolute(G_NN)
    
    diff = lagr_nn_difference(G_lgr, G_NN)
    print(f'G_NN shape is {G_NN.shape}')
    print(f'G_lgr shape is {G_lgr.shape}')
    print(f'l2 diff  shape is {diff.shape}')

    labels = ['Mode 1', 'Mode 10', 'Mode 20']
    save = f"G_diff_lagr_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_m{m}"
    ys = [diff[:, 0], diff[:, 9], diff[:, -1]]
    xlabel = 'Frequency (Hz)'
    ylabel = '$L_{2}$ normalised difference'
    plots.list_loglog_plot(freq_out, ys, xlabel, ylabel, labels, save)
    

def lagr_nn_difference(lgr, nn):
    diff_matrix = np.zeros(lgr.shape)
    lgr_norm = np.linalg.norm(lgr, axis=1)
    for j in range(lgr.shape[1]):
        diff_matrix[:, j] = (lgr[:, j] - nn[:, j])**2 / lgr_norm[j]
    return diff_matrix

if __name__ == '__main__':
    main()
