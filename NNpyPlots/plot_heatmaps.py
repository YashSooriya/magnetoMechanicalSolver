# import sys
# sys.path.insert(0 ,'../src_hpfem')
import funcs
import plots
import numpy as np
from scipy.io import loadmat


def save_data(solver, N_s, m, savename, customName):
    if len(customName) == 0:
        if m == 20:
            x_y = loadmat(f'arrays/X_y_{N_s}.mat')
        else:
            x_y = loadmat(f'arrays/X_y_{N_s}_m{m}.mat')
    else:
        x_y = loadmat(f'arrays/{customName}.mat')
    X = x_y['sample']
    y = x_y['y']
    print(f'X shape is {X.shape}.')
    print(f'y shape is {y.shape}.')
    losses, iterations, scores, layers, neurons = funcs.model_scores_layers_neurons(X, y, solver=solver, max_iter=100000)
    np.savez(f'arrays/heatmap_layers_neurons_{savename}', losses=losses,iterations=iterations,scores=scores,layers=layers,neurons=neurons)

load_data = lambda savename: np.load(f'arrays/heatmap_layers_neurons_{savename}.npz')


def plot_data(savename):
    data = load_data(savename)
    layers = data['layers']
    neurons = data['neurons']
    scores = data['scores']
    plots.plot_heatmap_layers_neurons(layers, neurons, scores, savename)
    
# -------------------------------------------------------------------------------------------------
load_xy = True
customLoad = False
solver = 'lbfgs'
n_s = 2324
m = 20
# -------------------------------------------------------------------------------------------------
if not customLoad:
    customName = ''
    if m == 20:
        savename = f'{solver}_s{n_s}'
    else:
        savename = f'{solver}_s{n_s}_m{m}'
elif customLoad:
    savename = input('Enter the x y data (without extension) in arrays/X_y_')
    customName = f"X_y_{savename}"
    savename = f'{solver}_{savename}'
print(customName)
if load_xy:
    save_data(solver, n_s, m, savename, customName)
plot_data(savename)
