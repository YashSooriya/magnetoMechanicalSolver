import plots
import numpy as np
from scipy.io import loadmat

m = 50
MULTIPLE = True

if MULTIPLE:
    save = f'singular_values_m{m}_logspace_comparison'
    Ns = [180, 90, 45, 23]
    singular_values_norm_list = []
    xs = []
    for N in Ns:
        data = loadmat(f'../data/SVDResult_Ns{N}_m{m}_5to5000Hz_logspace.mat')
        singular_values = np.diagonal(data['S'])
        singular_values_norm_list.append(singular_values /
                                         np.amax(singular_values))
        if N > 50:
            xs.append(np.arange(1, m + 1))
        else:
            xs.append(np.arange(1, N + 1))
    xticks = np.arange(0, xs[0].max()+1, 5)
    labels = ["$N_s = 180$", "$N_s = 90$", "$N_s = 45$", "$N_s = 23$"]
    plots.listxy_log_scatter_plot(xs, singular_values_norm_list, "Mode $m$",
                                  "Relative singular value", labels, save,
                                  xticks)
else:
    data = loadmat(f'../data/SVDResult_Ns{N}_m{m}_5to5000Hz_logspace.mat')        
    save = f'singular_values_m{m}_logspace'
    # data = loadmat(f'arrays/SVDResult_Ns180_m{m}.mat')
    # save = f'singular_values_m{m}'

    S = data['S']

    singular_values = np.diagonal(S)
    singular_values_norm = singular_values / np.amax(singular_values)

    print("singular values normalised: \n", singular_values_norm)

    singular_value_chosen = 20
    x = np.arange(1, singular_values.size + 1)
    xticks = np.arange(0, x.max()+1, 5)

    print(singular_values_norm[singular_value_chosen])
    plots.log_scatter_vline_plot(x, singular_values_norm, "Mode $m$", "Relative singular value", save, singular_value_chosen, xticks=xticks)







# data_tol10 = loadmat(f'arrays/SVDResult_Ns180_m{m}_tol1e-18.mat')
# data_tol18 = loadmat(f'arrays/SVDResult_Ns180_m{m}_tol1e-10.mat')
S = data['S']
# S_tol10 = data_tol10['S']
# S_tol18 = data_tol18['S']

singular_values = np.diagonal(S)
singular_values_norm = singular_values / np.amax(singular_values)
# singular_values_tol10 = np.diagonal(S_tol10)
# singular_values_norm_tol10 = singular_values_tol10 / np.amax(singular_values_tol10)
# singular_values_tol18 = np.diagonal(S_tol18)
# singular_values_norm_tol18 = singular_values_tol18 / np.amax(singular_values_tol18)

x = np.arange(1, singular_values.size + 1)
# y_list = [singular_values_norm_tol10, singular_values_norm, singular_values_norm_tol18]
labels = ['tolerance = $10^{-10}$','tolerance = $10^{-14}$', 'tolerance = $10^{-18}$']
xticks = np.arange(0, x.max()+1, 5)

# print("singular values normalised: \n", singular_values_norm)
# plots.list_log_scatter_plot(x, y, "Mode $m$", "Relative singular value", labels,"singular_values_tol",xticks=np.arange(0, x.max()+1, 5))
# singular_value_chosen = 20
# print(singular_values_norm[singular_value_chosen])
# plots.log_scatter_vline_plot(x, singular_values_norm, "Mode $m$", "Relative singular value", save, singular_value_chosen, xticks=xticks)
