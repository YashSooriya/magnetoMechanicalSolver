import numpy as np
import plots


N_s = 180

x = np.linspace(1, N_s, N_s)
freq_log = np.logspace(np.log10(5), np.log10(5000), N_s)
print(freq_log)


del_f1 = 10
del_f2 = 50

s1 = np.arange(10, 1000+1, del_f1)
s2 = np.arange(1000, 5000, del_f2)
freq_marcos = np.concatenate([s1, s2])

freqs = [freq_marcos, freq_log]
labels = ['Manually chosen snapshots', 'Log-space snapshots']
plots.list_scatter_plot(x, freqs, '#', 'Frequency', labels, 'marcos_logspace_freq', marker_size=3)

