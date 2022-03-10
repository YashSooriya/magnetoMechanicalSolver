import funcs
import plots
import numpy as np
from scipy.io import loadmat


# ns 359 180 90 45 23
# delf1 5 10 20 40 80
# delf2 25 50 100 200 400

# ---------------------------------------------------
del_f1 = 10
del_f2 = 50
m = 20
plot_mg = False
hist = True
imshow = False
# ---------------------------------------------------

# s1 = np.arange(10, 1000+1, del_f1)
# s2 = np.arange(1000, 5000, del_f2)
# samples = np.concatenate([s1, s2])

samples = np.logspace(np.log10(5),np.log10(5000),180)

# data = loadmat(f'arrays/SVDResult_Ns{samples.size}_m{m}.mat')
data = loadmat(f'../data/SVDResult_Ns{samples.size}_m{m}_5to5000Hz_logspace.mat')
G = data['G']

print(f'Snapshot frequency shapes is {samples.shape}')
print(f'G shape is {G.shape}')

G_amp = np.absolute(G)
save = f"s{samples.size}_m{m}_logged"
# save = f"s{samples.size}_m{m}"
if plot_mg:
    m_values = np.arange(1, m+1)
    ticks = np.arange(2, m+1, 2)
    plots.plot_m_G(m_values, G_amp, save, m, ticks)
else:
    plots.plot_freq_G(samples, G_amp, save, m)

if hist:
    plots.plot_G_hist(G_amp, save)
if imshow:
    plots.plot_G_imshow(G_amp, save)



