import funcs
import plots
import numpy as np
from scipy.io import loadmat
from loss import Loss


def main():
    Options = {'PODI': True, 'PODP': False, 'solver': 'trainbr', 'layers': 2, 'neurons': 16, 'm': 20,
               'per_mode': True, 'N_s': 180, 'network_no': 1}
    
    text_ops = {'text_PosX': 0.3, 'text_PosY': 0.25, 'text_diff': 0.08}


    MULTIPLE = True


    if MULTIPLE:
        modes = [1, 21, 39]

        Options['network_no'] = modes[0]
        loss_curve1 = Loss(Options)
        metrics,  save_name1 = loss_curve1.load()
        test_losses1 = metrics['test_losses']
        train_losses1 = metrics['train_losses']
        stop1 = metrics['stop'][0].replace("'","")
        epochs1 = metrics['epochs']
        
        Options['network_no'] = modes[1]        
        loss_curve2 = Loss(Options)
        metrics2, save_name2 = loss_curve2.load()
        test_losses2 = metrics2['test_losses']
        train_losses2 = metrics2['train_losses']
        stop2 = metrics2['stop'][0].replace("'","")
        epochs2 = metrics2['epochs']
        
        Options['network_no'] = modes[2]
        loss_curve3 = Loss(Options)
        metrics3, save_name3 = loss_curve3.load()
        test_losses3 = metrics3['test_losses']
        train_losses3 = metrics3['train_losses']
        stop3 = metrics3['stop'][0].replace("'","")
        epochs3 = metrics3['epochs']
        
        epochs = [epochs1, epochs2, epochs3]
        test_losses = [test_losses1, test_losses2, test_losses3]
        train_losses = [train_losses1, train_losses2, train_losses3]
        stop_reasons = [stop1, stop2, stop3]
        print(stop_reasons)
        plot_options = loss_curve1.plot_options()
        xlabel = plot_options['xlabel']
        ylabel = plot_options['ylabel']
        labels = ["mode 1", "mode 10", "mode 20"]
        save_name = save_name1 + f'_{modes[1]}_{modes[2]}'
        stop_labels = [x + ': ' + y for x,y in zip(labels, stop_reasons)]
        print(epochs)
        print(train_losses)

        save_name = save_name + f'_logspace'
        plots.plot_loss_curve(epochs, test_losses, train_losses, stop_labels, xlabel, ylabel, text_ops, save_name,
                              lgndlabels=labels,multiple=True)
        # plots.listxy_log_plot(epochs, train_losses, xlabel, ylabel, labels, save_name)

    else:
        loss_curve = Loss(Options)
        metrics, save_name = loss_curve.load()
        
        epochs = metrics['epochs']
        train_losses = metrics['train_losses']
        test_losses = metrics['test_losses']
        stop = metrics['stop'][0].replace("'","")
        best_train_loss = metrics['best_train_loss']
        best_epoch = metrics['best_epoch']
        
        plot_options = loss_curve.plot_options()
        xlabel = plot_options['xlabel']
        ylabel = plot_options['ylabel']
        labels = []

        plots.plot_loss_curve(epochs, test_losses, train_losses, stop, xlabel, ylabel, text_ops, save_name,
                              best_tr_loss=best_train_loss, best_tr_epoch = best_epoch)

if __name__ == '__main__':
    main()


