import numpy as np
import sys
sys.path.append('../src_hpfem')
sys.path.append('../NNPlots')
import plots
import funcs


def generate_simple_data(x):
    f1 = np.sin(x)
    f2 = np.cos(x)
    f3 = np.exp(x)
    f4 = x**2 + 3*x + 2
    return np.vstack((f1, f2, f3, f4)).T

funcslist = ['$\sin(x)$', '$\cos(x)$', '$\exp(x)$', '$x^2 + 3x + 2$']
activation = 'logistic'
layers = 1
neurons = 7
solver = 'lbfgs'
    
X = np.linspace(-1, 1, 100)
y = generate_simple_data(X)
print(y.shape)


#save_name = 'test_per_mode_funcs'
#plots.numpy_col_plot(X, y, 'x', '$f(x)$', funcslist, save_name)

predictX = np.linspace(-1, 1, 1000)
predictions = funcs.model_predict_per_mode(X, y, predictX, solver, activation, neurons, layers)
print(predictions.shape)

save_name_NN = f'test_per_mode_NN_{activation}_l{layers}_n{neurons}_funcs'
#plots.numpy_col_plot(predictX, predictions, 'x', '$f(x)$', funcslist_NN, save_name_NN)
plots.numpy_col_true_pred_plot(X, predictX, y, predictions, 'x', '$f(x)$', funcslist, save_name_NN)




