import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
import funcs
from normA import NormA
import seaborn as sns
from pandas import DataFrame
from matplotlib.colors import LogNorm
def main():
    Options = {'POD': 'PODI','neural_network': True, 'feed_forward_net': True, 'segmented': False, 'per_mode': True, 'N_s': 180, 'm': 30,
               'N_o': 40, 'x_axis_logged': False, 'No_segments': 20, 'layers': 1, 'chosen_sample': 'log space',
               'neurons': 1, 'solver': 'adam', 'activation': 'identity', 'ffn_layers': 1, 'ffn_neurons': 16, 'ffn_solver': 'trainbr',
               'custom_saveload': False, 'save_name': '', 'load_name': ''}
    MIN_NEURONS = 8
    MAX_NEURONS = 32
    MAX_LAYERS = 5
    neurons_arr = np.logspace(int(np.log2(MIN_NEURONS)), int(np.log2(MAX_NEURONS)),num=int(np.log2(MAX_NEURONS))-int(np.log2(MIN_NEURONS))+1, base=2, dtype=int)
    layers_arr = np.arange(1, MAX_LAYERS+1, 1)
    # freq_out = np.linspace(15,15+(Options['N_o']-1)*10, Options['N_o'])
    # freq_out = np.linspace(5, 105, N_o)
    freq_out = np.linspace(15, 5000, Options['N_o'])
    Options['freqout'] = freq_out


    y_lagr = load_lagrange(Options)

    print(layers_arr)
    print(neurons_arr)
    
    scores = np.zeros((len(layers_arr), len(neurons_arr)))
    for i, layers in enumerate(layers_arr):
        Options['ffn_layers'] = layers
        for j, neurons in enumerate(neurons_arr):
            Options['ffn_neurons'] = neurons
            normA_nn = NormA(Options)
            y_nn, _ = normA_nn.load()
            difference = lagr_nn_difference(y_lagr, y_nn, average=True)
            scores[i, j] = difference

    # scores = normalise_data(scores)
    # scores = 1- scores
    scores = 1/scores
    print(scores)
    savename = f'l{layers_arr[0]}tol{layers_arr[-1]}_n{neurons_arr[0]}ton{neurons_arr[-1]}_logspace'
    # plots.plot_heatmap_layers_neurons_ylog(layers_arr, neurons_arr, scores, savename)

    df = DataFrame(scores.T, index=neurons_arr,columns=layers_arr)
    sns.heatmap(df, cmap=plt.cm.BuPu,cbar=True, norm=LogNorm(),cbar_kws={'label':'Score'})
    plt.xlabel("No. of layers")
    plt.ylabel("No. of neurons")
    # plt.yscale('log',basey=2)
    plt.gca().invert_yaxis()
    plt.savefig(f'figures/heatmap_layers_neurons_{savename}.pdf')
    plt.savefig(f'figures/heatmap_layers_neurons_{savename}.eps')
    print(f'saved to figures/heatmap_layers_neurons_{savename}')
    plt.show()
    


def load_lagrange(Options):
    old_m = Options['m']
    Options['neural_network'] = False
    Options['m'] = 20
    normA_lagr = NormA(Options)
    y_lagr, _ = normA_lagr.load()
    Options['neural_network'] = True
    Options['m'] = old_m    
    return y_lagr

    
def normalise_data(diff):
    return (diff - np.nanmin(diff)) / (np.nanmax(diff) - np.nanmin(diff))
        
    
def lagr_nn_difference(y1, y2, average=True):
    norm_diff_shields = np.zeros((3, y1[0].shape[0]))
    for i in range(3):
        diff = np.square(y1[i] - y2[i]) / np.nanmax(y1[i])
        norm_diff_shields[i, :]= diff.ravel()
            

    if average:
        ret = np.nanmean(norm_diff_shields, axis=None)
    else:
        ret = norm_diff_shields
    return ret


if __name__ == '__main__':
    main()
