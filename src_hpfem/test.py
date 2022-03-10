import numpy as np
import os
import funcs
import matplotlib.pyplot as plt

# plot line colours, 10 colors all together
CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'
CB91_Red = '#EA1200'
CB91_Navy = '#1F618D'
CB91_Grey = '#797D7f'
CB91_Black = '#17202A'


color_list = [CB91_Blue, CB91_Pink, CB91_Green, CB91_Amber,
              CB91_Purple, CB91_Violet, CB91_Red, CB91_Navy, CB91_Grey, CB91_Black]

# X = np.linspace(15, 5000, 20)
# del_f1 = 10
# del_f2 = 50
# X1 = np.arange(10, 1000+1, del_f1)
# X2 = np.arange(1000, 5000, del_f2)
# X = np.concatenate((X1, X2))
# print(X.shape)
# y = np.random.randint(1, high=50, size=(20, 50))
#y = np.random.rand(X.shape[0], 50)
# print(y)


# simple functions
def generate_simple_data(x):
    f1 = np.sin(x)
    f2 = np.cos(x)
    f3 = np.exp(x)
    f4 = x**2 + 3*x + 2
    return np.vstack((f1, f2, f3, f4)).T

def plot_data(x, y, xlabel, ylabel):
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)
    plt.rcParams['lines.linewidth'] = 2
    funcslist = ['$\sin(x)$', '$\cos(x)$', '$\exp(x)$', '$x^2 + 3x + 2$']
    for index in range(y.shape[1]):
        plt.plot(x, y[:, index], label=f'{funcslist[index]}')
    plt.legend()
    plt.show()
    
X = np.linspace(-1, 1, 100)
y = generate_simple_data(X)
print(y.shape)
plot_data(X, y, 'x', '$f(x)$')
# print(y)
predictX = np.linspace(-1, -1, 1000)
predictions = funcs.model_predict_per_mode(X, y, predictX, 'adam', 'identity', 1, 1)
print(predictions)
print(predictions.shape)
