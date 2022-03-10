import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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

plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)
plt.rcParams['lines.linewidth'] = 2
# plt.rcParams['text.usetex'] = True

def plot_loss_curve(epochs, test_loss, train_loss, stop_reasons, xlabel, ylabel, text_ops, save_name,
                    lgndlabels=[], best_tr_loss=0, best_tr_epoch=0, multiple=False):
    plt.figure()
    text_posX = text_ops['text_PosX']
    text_posY = text_ops['text_PosY']
    text_diff = text_ops['text_diff']
    if multiple:
        assert(len(test_loss) == len(train_loss))
        n_plots = len(test_loss)
        for i in range(n_plots):
            # print(f'{epochs[i]=}')
            # print(f'{test_loss[i]=}')
            # print(f'{train_loss[i]=}')
            plt.plot(epochs[i], train_loss[i], label="Train: {}".format(lgndlabels[i]))
            plt.plot(epochs[i], test_loss[i], label="Test: {}".format(lgndlabels[i]))
            plt.text(text_posX, text_posY - text_diff*i , stop_reasons[i], bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'red',
                                                            'pad':5}, ha='center', va='center', transform=plt.gca().transAxes)
    else:
        plt.plot(epochs, train_loss, label='Train')
        plt.plot(epochs, test_loss, label='Test')
        plt.axhline(y=best_tr_loss, color='k', linestyle='--', label=f'Best train = {best_tr_loss[0,0]:.5}', linewidth=1)
        plt.axvline(x=best_tr_epoch, color='k', linestyle='--', linewidth=1)
#        plt.Circle((best_tr_epoch, best_tr_loss), 0.0001, color='blue')

        plt.text(text_posX, text_posY - text_diff, stop_reasons, bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'red',
                                                                             'pad':5}, ha='center', va='center', transform=plt.gca().transAxes)

    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=0)
    plt.yscale("log")
    plt.savefig(f"figures/{save_name}.pdf")
    print(f'Saved figure to figures/{save_name}.pdf')
    plt.show()
        
def numpy_col_true_pred_plot(x_true, x_pred, true_array, pred_array, xlabel, ylabel, lgndlabels, save_name):
    ncols = true_array.shape[1]
    for i in range(ncols):
        plt.figure()
        plt.plot(x_true, true_array[:, i], label=f'True {lgndlabels[i]}')
        plt.plot(x_pred, pred_array[:, i], label=f'Predicted {lgndlabels[i]}')
        plt.legend()
        plt.savefig(f'figures/{save_name}_{lgndlabels[i]}.pdf')
        print(f'figures/{save_name}_{lgndlabels[i]}.pdf')
        plt.show()

def list_scatter_line_log_plot(x1,x2, ylist, xlabel, ylabel, lgndlabels, save, marker_size=6):
    plt.figure()
    for i in range(len(ylist)):
        plt.plot(x1, ylist[i][0], linestyle="None", marker="x", label=lgndlabels[i][0], markersize=marker_size)
        plt.plot(x2, ylist[i][1], label=lgndlabels[i][1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.yscale("log")
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()

def list_scatter_line_loglog_plot(x1,x2, ylist, xlabel, ylabel, lgndlabels, save, marker_size=6):
    plt.figure()
    for i in range(len(ylist)):
        plt.plot(x1, ylist[i][0], linestyle="None", marker="x", label=lgndlabels[i][0], markersize=marker_size)
        plt.plot(x2, ylist[i][1], label=lgndlabels[i][1])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()


def list_scatter_line_log_multiple_plot(x1,x2,ylist,xlabel, ylabel, lgndlabels, save, marker_size=6):
    plt.figure()
    plt.plot(x2, ylist[0], linestyle="None", marker="x", label=lgndlabels[0], markersize=marker_size)
    for i in range(1, len(ylist)):
        plt.plot(x1, ylist[i], label = lgndlabels[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.yscale("log")
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()



        
def numpy_col_plot(x, y, xlabel, ylabel, lgndlabels, save):
    plt.figure()
    for index in range(y.shape[1]):
        plt.plot(x, y[:, index], label=lgndlabels[index])
    plt.legend()
    plt.savefig(f'figures/{save}.pdf')
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()


def simple_log_plot(x, y, xlabel, ylabel, save, xticks=[]):
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()

def list_log_diff_plot(x, ylist, xlabel, ylabel, lgndlabels, save, diff,xticks=[]):
    plt.figure(figsize=(12,9), dpi=100)
    for i in range(len(ylist)):
        plt.plot(x, ylist[i], label = lgndlabels[i])
    plt.plot(x, diff, 'k--', label = 'Average Difference')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()

def listxy_log_plot(xlist, ylist, xlabel, ylabel ,lgndlabels, save, xticks=[], lgndlocation=0):
    for i in range(len(ylist)):
        plt.plot(xlist[i], ylist[i], label = lgndlabels[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=lgndlocation)
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()

def listxy_loglog_plot(xlist, ylist, xlabel, ylabel ,lgndlabels, save, xticks=[], lgndlocation=0):
    for i in range(len(ylist)):
        plt.plot(xlist[i], ylist[i], label = lgndlabels[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=lgndlocation)
    plt.yscale("log")
    plt.xscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()


    
def list_log_plot(x, ylist, xlabel, ylabel, lgndlabels, save, xticks=[], lgndlocation=0):
    for i in range(len(ylist)):
        plt.plot(x, ylist[i], label = lgndlabels[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=lgndlocation)
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()

def list_loglog_plot(x, ylist, xlabel, ylabel, lgndlabels, save, xticks=[], lgndlocation=0):
    for i in range(len(ylist)):
        plt.plot(x, ylist[i], label = lgndlabels[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(loc=lgndlocation)
    plt.yscale("log")
    plt.xscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()



def list_xylog_plot(x, ylist, xlabel, ylabel, lgndlabels, save, xticks=[]):
    for i in range(len(ylist)):
        plt.plot(x, ylist[i], label = lgndlabels[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.yscale("log")
    plt.xscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()



def log_scatter_plot(x, y, xlabel, ylabel, save, xticks=[]):
    plt.plot(x, y, linestyle="None", marker="x")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()

def log_scatter_vline_plot(x, y, xlabel, ylabel, save, vline,xticks=[]):
    plt.plot(x, y, linestyle="None", marker="x")
    plt.vlines(vline, 0, np.max(y), 'r', 'dashed')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f'Saved figure to figures/{save}.pdf')
    plt.show()

    
def list_log_scatter_plot(x, ylist, xlabel, ylabel, lgndlabels, save, xticks=[]):
    for i in range(len(ylist)):
        plt.plot(x, ylist[i], label = lgndlabels[i], linestyle="None", marker="x")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f"Saved to figures/{save}.pdf")
    plt.show()

def listxy_log_scatter_plot(x, ylist, xlabel, ylabel, lgndlabels, save, xticks):
    for i in range(len(ylist)):
        plt.plot(x[i], ylist[i], label = lgndlabels[i], linestyle="None", marker="x")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.yscale("log")
    if len(xticks) != 0:
        plt.xticks(xticks)
    plt.savefig(f"figures/{save}.pdf")
    print(f"Saved to figures/{save}.pdf")
    plt.show()
    
def simple_plot(x, y, xlabel, ylabel, save):
    plt.plot(x, y)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(f"figures/{save}.pdf")
    print(f"Saved to figures/{save}.pdf")
    plt.show()

def simple_scatter_plot(x, y, xlabel, ylabel, save):
    plt.plot(x, y, linestyle="None", marker="x")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(f"figures/{save}.pdf")
    print(f"Saved to figures/{save}.pdf")
    plt.show()


    
def list_plot(x, ylist, xlabel, ylabel, lgndlabels, save):
    for i in range(len(ylist)):
        plt.plot(x, ylist[i], label = lgndlabels[i])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(f"figures/{save}.pdf")
    plt.show()

def list_scatter_plot(x, ylist, xlabel, ylabel, lgndlabels, save, marker_size=6):
    for i in range(len(ylist)):
        plt.plot(x, ylist[i], label = lgndlabels[i], linestyle="None", marker="x", markersize=marker_size)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.savefig(f"figures/{save}.pdf")
    print(f"Saved to figures/{save}.pdf")

    plt.show()

def plot_freq_G(freq, G, savename, m):
    max_lines = 10
    f1 = plt.figure()
    f2 = plt.figure()
    ax1 = f1.add_subplot(111)
    ax2 = f2.add_subplot(111)
    lines = np.arange(G.shape[1])

    for j in lines:
        if j < max_lines:
            ax1.plot(freq, G[:,j], label=f"$\mathbf{{G}}_{{i{j}}}$")
        else:
            ax2.plot(freq, G[:,j], label=f"$\mathbf{{G}}_{{i{j}}}$")
            
    ax1.set_yscale("log")
    ax1.legend()
    ax1.set_xlabel("Frequency $i$ (Hz)")
    ax1.set_ylabel("|$\mathbf{G_{ij}}$|")
    ax1.set_xscale("log")
    f1.savefig(f"figures/G_lineplot_1-10_{savename}.pdf")
    if m > max_lines:
        ax2.set_yscale("log")
        ax2.legend()
        ax2.set_xlabel("Frequency $i$ (Hz)")
        ax2.set_ylabel("|$\mathbf{G_{ij}}$|")
        ax2.set_xscale("log")
        f2.savefig(f"figures/G_lineplot_11-20_{savename}.pdf")
    plt.show()
        
def plot_m_G(m, G, savename, n_freqs, xticks):
    max_lines = 10
    f1 = plt.figure()
    f2 = plt.figure()
    ax1 = f1.add_subplot(111)
    ax2 = f2.add_subplot(111)
    lines = np.arange(G.shape[0])
    markers = ["x", "*", "^", "v", ">", "<", "1", "2", "3", "4", "+", "D", "."]

    r = 0
    q = 0
    # or in lines
    for i in range(20):
        if i < max_lines:
            ax1.plot(m, G[i,:], label=f"$\mathbf{{G}}_{{{i}j}}$", linestyle="None", marker=markers[r])
            r += 1
        else:
            ax2.plot(m, G[i,:], label=f"$\mathbf{{G}}_{{{i}j}}$", linestyle="None", marker=markers[q])
            q += 1
            
    ax1.set_yscale("log")
    ax1.legend()
    ax1.set_xlabel("Mode $j$")
    ax1.set_ylabel("|$\mathbf{G_{ij}}$|")
    ax1.set_xticks(xticks)
    f1.savefig(f"figures/G_lineplot_xm_1-10_{savename}.pdf")
    if n_freqs > max_lines:
        ax2.set_yscale("log")
        ax2.legend()
        ax2.set_xlabel("Mode $j$")
        ax2.set_ylabel("|$\mathbf{G_{ij}}$|")
        ax2.set_xticks(xticks)
        f2.savefig(f"figures/G_lineplot_xm_11-20_{savename}.pdf")
    plt.show()

def plot_G_hist(G, savename):
    plt.figure()
    logbins = np.logspace(np.log10(G.min()), np.log10(G.max()), 10)
    bins = np.linspace(G.min(), G.max(), 10)
    print(logbins)
    G = G.flatten()
    plt.hist(G, bins=logbins, rwidth=0.7, alpha = 0.7, )
    plt.xscale("log")
    # plt.xticks(logbins)
    plt.xlabel("$|\mathbf{G}_{ij}|$")
    plt.ylabel("Frequency")
    # plt.xticks(fontsize=15)
    # plt.yticks(fontsize=15)
    plt.savefig(f"figures/G_hist_{savename}.pdf")
    print(f"Saved to figures/G_hist_{savename}.pdf")
    plt.show()

def plot_G_imshow(G, savename):
    plt.figure()
    plt.imshow(G, interpolation='none', cmap='BuPu', norm=LogNorm(vmin=1e-6, vmax=1))
    plt.colorbar(label='$|\mathbf{G}_{ij}|$')
    plt.xlabel("Modes $j$")
    plt.ylabel("Frequencies $i$")
    # plt.xticks(fontsize=15)
    # plt.yticks(fontsize=15)
    plt.savefig(f"figures/G_imshow_{savename}.pdf")
    print(f"Saved to figures/G_imshow_{savename}.pdf")
    plt.show()



# plot loss agianst no. of iterations for varying number of neurons
def plot_neurons_loss(neuron_iterations, loss, neurons, savename):
    i = 0
    for iterations, loss in zip(neuron_iterations.T, loss.T):
        plt.plot(iterations, loss, label="hidden neurons = {}".format(neurons[i]))
        i += 1
    plt.legend()
    plt.xlabel("No. of iterations")
    plt.ylabel("Relative loss")
    plt.savefig(f"figures/neurons_loss_{savename}.pdf")
    plt.savefig(f"figures/neurons_loss_{savename}.eps")
    plt.show()

def plot_neurons_log_loss(neuron_iterations, loss, neurons, savename):
    plt.figure()
    i = 0
    for iterations, loss in zip(neuron_iterations.T, loss.T):
        plt.plot(iterations, loss, label="hidden neurons = {}".format(neurons[i]))
        i += 1
    plt.legend()
    plt.xlabel("No. of iterations")
    plt.ylabel("Loss")
    plt.yscale('log')
    plt.savefig(f"figures/neurons_logloss_{savename}.pdf")
    plt.savefig(f"figures/neurons_logloss_{savename}.eps")
    plt.show()
    
def plot_neurons_score(neurons, score, savename):
    plt.figure()
    plt.rcParams['axes.prop_cycle'] = plt.cycler(color=color_list)
    plt.rcParams['lines.linewidth'] = 2
    plt.plot(neurons, score)
    plt.xlabel("No. of neurons")
    plt.ylabel("$R^2$ score")
    plt.savefig(f"figures/neurons_score_{savename}.pdf")
    plt.savefig(f"figures/neurons_score_{savename}.eps")
    plt.show()

def plot_heatmap_layers_neurons(layers, neurons, scores_array, savename):
    # Nx = len(layers)
    # Ny = len(neurons)
    # x = np.linspace(0, max(layers)-1, Nx)
    # y = np.linspace(0, max(neurons)-1, Ny)
    scores_array = scores_array.T
    x = np.array(layers)
    y = np.array(neurons)
    [xx, yy] = np.meshgrid(x, y)
    print(f'pcolor xx shape: {xx.shape}')
    print(f'pcolor yy shape: {yy.shape}')
    plt.pcolor(xx, yy, scores_array, cmap='BuPu', shading='auto', vmin=-0.1)
#    print(np.where(scores_array.T == max(scores_array.T)))
    max_indicies = np.unravel_index(np.argmax(scores_array, axis=None), scores_array.shape)
    # print(max_indicies)
    print(f'layer with heighest score: {x[max_indicies[1]]}')
    print(f'neuron with heighest score: {y[max_indicies[0]]}')
    plt.xlabel("No. of layers")
    plt.ylabel("No. of neurons")
    plt.colorbar(label='$R^2$ score', ticks=np.arange(-0.1, np.amax(scores_array), 0.05))
    plt.savefig(f'figures/heatmap_layers_neurons_{savename}.pdf')
    plt.savefig(f'figures/heatmap_layers_neurons_{savename}.eps')
    plt.show()

def plot_heatmap_layers_neurons_ylog(layers, neurons, scores_array, savename):
    # Nx = len(layers)
    # Ny = len(neurons)
    # x = np.linspace(0, max(layers)-1, Nx)
    # y = np.linspace(0, max(neurons)-1, Ny)
    scores_array = scores_array.T
    x = np.array(layers)
    y = np.array(neurons)
    [xx, yy] = np.meshgrid(x, y)
    print(f'pcolor xx shape: {xx.shape}')
    print(f'pcolor yy shape: {yy.shape}')
    plt.pcolor(xx, yy, scores_array, cmap='BuPu', shading='auto')
#    print(np.where(scores_array.T == max(scores_array.T)))
    max_indicies = np.unravel_index(np.argmax(scores_array, axis=None), scores_array.shape)
    # print(max_indicies)
    print(f'layer with heighest score: {x[max_indicies[1]]}')
    print(f'neuron with heighest score: {y[max_indicies[0]]}')
    plt.xlabel("No. of layers")
    plt.ylabel("No. of neurons")
    plt.yticks(neurons)
    plt.yscale('log',basey=2)
#    plt.colorbar(label='Score', ticks=np.arange(0, np.amax(scores_array), 0.05))
    plt.colorbar(label='Score')
    plt.savefig(f'figures/heatmap_layers_neurons_{savename}.pdf')
    plt.savefig(f'figures/heatmap_layers_neurons_{savename}.eps')
    plt.show()
    
    
    
