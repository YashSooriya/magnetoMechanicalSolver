import numpy as np
from scipy.linalg.decomp_svd import null_space
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from scipy.io import savemat

def test(A):
    B = np.array([1,2,3])
    return A + B

# prediction using a NN model with a given set of hyperparameters
def model_predict(X, y,predict_data, solver, activation, neurons_in_each_layer, layers, save_split):
    # since there is a single feature in X
    X = X.reshape(-1, 1)
    predict_data = predict_data.reshape(-1, 1)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
    # print(f'X_train shape is {X_train.shape}')
    
    if save_split:
        savedic = {"X_train": X_train, "X_test": X_test, "y_train": y_train, "y_test": y_test}
        savemat("data/NN_train_test_split.mat", savedic)
        
    X_train, X_test, y_train, y_test, scalerX, scalerY = scale_train_test(X_train, X_test, y_train, y_test)
    predict_data_s = scalerX.transform(predict_data)

    kwargs = {"solver": solver, "activation": activation, 
    "hidden_layer_sizes": tuple([neurons_in_each_layer]*layers), 
            "random_state": 42, "max_iter": 10000}

    regr = MLPRegressor(**kwargs)
    regr = regr.fit(X_train, y_train)
    score = regr.score(X_test, y_test)

    prediction_s = regr.predict(predict_data_s)
    prediction = scalerY.inverse_transform(prediction_s)
    return prediction, score, X_train

# prediction using a NN model for each segment of X with a given set of hyperparameters
def model_predict_segments(X, y,predict_data, solver, activation, neurons_in_each_layer, layers, no_segments):
    # splitting into segments

    X_segments = np.array_split(X, no_segments)
    y_segments = np.array_split(y, no_segments)
    xdiff = X[1] - X[0]

    # X_segments = np.split(X, X.size/segment_size)
    # y_segments = np.split(y, y.shape[0]/segment_size)
    min_segments = np.zeros(no_segments)
    max_segments = np.zeros(no_segments)
    networks = []
    scalerXs = []
    scalerYs = []
#    scores = np.zeros(len(X_segments))
    i = 0
    for X_segment, y_segment in zip(X_segments, y_segments):
        
        # since there is a single feature in X
        X_segment = X_segment.reshape(-1, 1)
        predict_data = predict_data.reshape(-1, 1)

        # All data is training
        X_train = X_segment
        y_train = y_segment
        # print(f'X_train shape is {X_train.shape}')
        
        scalerX = StandardScaler().fit(X_train.reshape(-1, 1))
        X_train = scalerX.transform(X_train)

        scalerY = StandardScaler().fit(y_train)
        y_train = scalerY.transform(y_train)

        scalerXs.append(scalerX)
        scalerYs.append(scalerY)
        kwargs = {"solver": solver, "activation": activation, 
                "hidden_layer_sizes": tuple([neurons_in_each_layer]*layers), 
                "random_state": 42, "max_iter": 10000}
    
        regr = MLPRegressor(**kwargs)
        regr = regr.fit(X_train, y_train)
        # score = regr.score(X_test, y_test)
        # scores[i] = score
        networks.append(regr)
        min_segments[i] = max_segments[i-1]
        max_segments[i] = np.max(X_segment)
        i += 1
        

    # min_segments = np.min(X_segments, axis=1) - xdiff
    # max_segments = np.max(X_segments, axis=1)
#    print(min_segments)
#    print(max_segments)
    predictions = np.zeros((predict_data.shape[0], y.shape[1]))
    i = 0
    # predict_data = predict_data.reshape(-1, 1)
    for data_point in predict_data:
        if data_point <= min_segments[0]:
            index = 0
        elif data_point >= max_segments[max_segments.size-1]:
            index = max_segments.size - 1
        else:
            msk = (data_point <= max_segments) & (data_point >= min_segments)
 #           print(msk)
            index = np.argwhere(msk)[0, 0]
        
        scalerX = scalerXs[index]
        scalerY = scalerYs[index]
        regr = networks[index]

        data_point = data_point.reshape(-1, 1)
        predict_data_s = scalerX.transform(data_point)
        prediction_s = regr.predict(predict_data_s)
        prediction = scalerY.inverse_transform(prediction_s)
        predictions[i, :] = prediction
        i += 1
    return predictions

def model_predict_per_mode(X, y,predict_data, solver, activation, neurons_in_each_layer, layers):
    n_modes = y.shape[1]
    X = X.reshape(-1, 1)
    
    # All data is training
    X_train = X
    scalerX = StandardScaler().fit(X_train.reshape(-1, 1))
    X_train = scalerX.transform(X_train)

    predictions = np.zeros((predict_data.shape[0], n_modes))
    
    predict_data = predict_data.reshape(-1, 1)
    predict_data = scalerX.transform(predict_data)

    for mode in range(n_modes):
        y_train = y[:, mode]
        y_train = y_train.reshape(-1, 1)
        scalerY = StandardScaler().fit(y_train)
        y_train = scalerY.transform(y_train)
        y_train = y_train.reshape(-1)

        kwargs = {"solver": solver, "activation": activation, 
                "hidden_layer_sizes": tuple([neurons_in_each_layer]*layers), 
                "random_state": 42, "max_iter": 10000}
    
        regr = MLPRegressor(**kwargs)
        regr = regr.fit(X_train, y_train)
        # best_loss = regr.best_loss_
        
        prediction = regr.predict(predict_data)
        prediction = prediction.reshape(-1, 1)
        prediction = scalerY.inverse_transform(prediction)
        prediction = prediction.reshape(-1)
        predictions[:, mode] = prediction

    return predictions


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

print("Successfully loaded funcs.py")
