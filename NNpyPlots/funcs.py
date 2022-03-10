import numpy as np
from scipy.linalg.decomp_svd import null_space
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from scipy.io import savemat

# loop over from min neurons to max neurons in steps and get the scores for each neuron size, for a given layer number
def model_scores_neuron(X, y,solver='adam',activation='relu', min_neurons=5,max_neurons=10, step=2,layers=1, max_iter = 10000):
    # since there is a single feature in X
    X = X.reshape(-1, 1)
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    X_train, X_test, y_train, y_test, scalerX, scalerY = scale_train_test(X_train, X_test, y_train, y_test)

    neurons = np.arange(min_neurons, max_neurons+1, step)
    relative_loss_denom = np.mean(np.square(y_train))
    losses_arr = np.zeros(neurons.size,dtype=object)
    iterations_arr = np.zeros(neurons.size,dtype=object)
    scores_arr = np.zeros(neurons.size)

    i = 0
    for neuron in neurons:
        kwargs = {"solver": solver, "activation": activation, 
        "hidden_layer_sizes": tuple([neuron]*layers), 
                  "random_state": 42, "max_iter": max_iter}
        regr = MLPRegressor(**kwargs)
        regr = regr.fit(X_train, y_train)
        score = regr.score(X_test, y_test)
        
        if solver == 'adam' or solver == 'sgd':
            n_iterations = regr.n_iter_
            loss = regr.loss_curve_
            relative_loss = loss / relative_loss_denom
            iterations = np.linspace(1,n_iterations, len(loss))
            losses_arr[i] = relative_loss
            iterations_arr[i] = iterations
        
        scores_arr[i] = score
        i += 1 
    return losses_arr, iterations_arr, scores_arr, neurons


# loop over both neurons and layers, the neuron number is the same in every layer.
def model_scores_layers_neurons(X, y, solver='adam', activation='relu', min_neurons=1, max_neurons=30, neuron_step=1, min_layers=1, max_layers=10, layers_step=1, max_iter=10000):
    # reshape due to single feature
    X = X.reshape(-1, 1)

    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    X_train, X_test, y_train, y_test, scalerX, scalerY = scale_train_test(X_train, X_test, y_train, y_test)

    neurons = np.arange(min_neurons, max_neurons+1, neuron_step)
    layers = np.arange(min_layers, max_layers+1, layers_step)
    
    losses_arr = np.zeros((max_layers, max_neurons),dtype=object)
    iterations_arr = np.zeros((max_layers, max_neurons))
    scores_arr = np.zeros((max_layers, max_neurons))

    i = 0
    for layer in layers:
        j = 0
        for neuron in neurons:
            kwargs = {"solver": solver, "activation": activation, 
            "hidden_layer_sizes": tuple([neuron]*layer), 
                      "random_state": 42, "max_iter": max_iter}
            regr = MLPRegressor(**kwargs)
            regr = regr.fit(X_train, y_train)
            
            if solver == 'adam' or solver == 'sgd':
                n_iterations = regr.n_iter_
                iterations_arr[i, j] = n_iterations
                loss = regr.loss_curve_
                losses_arr[i, j] = loss
                
            score = regr.score(X_test, y_test)
            scores_arr[i, j] = score
            j += 1
        i += 1
        print(f'layer {layer} complete.')
    return losses_arr, iterations_arr, scores_arr, layers, neurons

# scale both the X and y train and test data to single variance and 0 mean    
def scale_train_test(X_train, X_test, y_train, y_test):
    scalerX = StandardScaler().fit(X_train.reshape(-1, 1))
    scaled_X_train, scaled_X_test = scale(X_train, X_test, scalerX)

    scalerY = StandardScaler().fit(y_train)
    scaled_y_train, scaled_y_test = scale(y_train, y_test, scalerY)
    
    return scaled_X_train, scaled_X_test, scaled_y_train, scaled_y_test, scalerX, scalerY
        
def scale(X1, X2, scaler, unscale=False):
    if unscale:
        X1_s = scaler.inverse_transform(X1)
        X2_s = scaler.inverse_transform(X2)
    else:
        X1_s = scaler.transform(X1)
        X2_s = scaler.transform(X2)
    return X1_s, X2_s
