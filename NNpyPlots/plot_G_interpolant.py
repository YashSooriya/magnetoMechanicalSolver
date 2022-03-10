import plots
import numpy as np
from scipy.io import loadmat


def main():
    Options = {'N_s': 23, 'm': 20, 'N_o': 40, 'chosen_sample': 'log space',
               'neurons': 1, 'solver': 'adam', 'activation': 'identity',
               'ffn_layers': 3, 'ffn_neurons': 16, 'ffn_solver': 'trainbr',
               'min_grad': 1e-10, 'reg': 0}

    globals().update(Options)

    samples = np.logspace(np.log10(5), np.log10(5000), N_s)
    data_G = loadmat(f'../data/SVDResult_Ns{N_s}_m{m}_5to5000Hz_logspace.mat')
    G = data_G['G']
    data_int = loadmat(f'data/GCurve/GCurve_{ffn_solver}_l{ffn_layers}_n{ffn_neurons}_m{m}_Ns{N_s}_reg{reg}.mat')
    freq_out = np.linspace(15, 5000, N_o)
    G_int = data_int['new_G']

    print(f'G shape is {G.shape}')
    print(f'G_int shape is {G_int.shape}')

    G = np.absolute(G)
    G_int = np.absolute(G_int)

    x = [samples, freq_out, samples, freq_out, samples, freq_out]
    y = [G[:, 0], G_int[:, 0], G[:, 9], G_int[:, 9], G[:, -1], G_int[:, -1]]
    
    ys = [[G[:, i], G_int[:, i]] for i in [0, 9, -1]]
    
    labels = ['SVD $\mathbf{G}_1$', 'Interpolated $\mathbf{G}_1$', 'SVD $\mathbf{G}_{10}$', 'Interpolated $\mathbf{G}_{10}$',
              'SVD $\mathbf{G}_{20}$', 'Interpolated $\mathbf{G}_{20}$']
    labels2 = [[f'SVD $\mathbf{{G}}_{{{i}}}$', f'Interpolated $\mathbf{{G}}_{{{i}}}$'] for i in [0, 9,20]]
    print(labels2)
    save = f"G_int_comparison_l{ffn_layers}_n{ffn_neurons}_Ns{N_s}_m{m}"

    xlabel = 'Frequency (Hz)'
    ylabel = '$||\mathbf{{G}}_{{j}}||$'
    # plots.listxy_loglog_plot(x, y, xlabel, ylabel, labels, save)
    plots.list_scatter_line_loglog_plot(samples, freq_out, ys, xlabel, ylabel, labels2, save, marker_size=4)
    






if __name__ == '__main__':
    main()
