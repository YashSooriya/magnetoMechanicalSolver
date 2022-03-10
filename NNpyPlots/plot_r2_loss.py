# import sys
# sys.path.insert(0 ,'../src_hpfem')
import funcs
import plots
import numpy as np
from scipy.io import loadmat

x_y = loadmat('arrays/X_y.mat')
X = x_y['sample']
y = x_y['G_complex_real']
print(f'X shape is {X.shape}.')
print(f'y shape is {y.shape}.')

loss, iterations, r2scores, neurons = funcs.model_scores_neuron(X, y, solver='adam', activation='relu',max_neurons=30, step=5)
plots.plot_neurons_loss(iterations, loss, neurons, 'adam')
plots.plot_neurons_score(neurons, r2scores, 'adam')


# loss, iterations, r2scores, neurons = funcs.model_scores_neuron(X, y, solver='lbfgs', activation='relu',max_neurons=30, step=5, max_iter=20000)
# plots.plot_neurons_score(neurons, r2scores, 'lbfgs')





# plots.plot_neurons_log_loss(iterations, loss, neurons)



