function [new_q] = neuralNetworkSVDPODP(freqSample, H, S, G, trainSet, freqOut, Options)
    % D = H * S * G';
    % trainSet = H'* D;
    % fprintf('Size of D is (%d, %d) \n', size(D))
    fprintf('Size of trainSet is (%d, %d) \n', size(trainSet))
    % trainSet = (inv(S)*trainSet)';
    trainSet = trainSet';
    % (G - trainSet)

    if Options.feedForwardNet==1
        y = transpose(complex2cols(conj(trainSet)));
        % y = (y - mn)/st;
        % y(40, :)
        fprintf('Size of y is (%d, %d) \n', size(y))
        noNetworks = size(y, 1);
        solver = Options.ffnSolver;
        hiddenLayers = Options.ffnLayers;
        hiddenNeurons = Options.ffnNeurons;
        predictions = zeros(noNetworks, size(freqOut, 2));
        hiddenSizes = ones(1, hiddenLayers)*hiddenNeurons;
        
        if Options.ffnPerMode==1
            saveFolder = sprintf('NNdata/loss_curves/PODP/loss_curves_%s_l%d_n%d_m%d_Ns%d/', solver, hiddenLayers, hiddenNeurons, noNetworks/2, size(trainSet, 2));
            if ~exist(saveFolder, 'dir')
                mkdir(saveFolder)
            end
            startTrain = tic;
            for i = 1:noNetworks
                if Options.functionFitNet == 1
                    net = fitnet(hiddenSizes, solver);
                else
                    net = feedforwardnet(hiddenSizes, solver);
                end
                
                % Setup Division of Data for Training, Validation, Testing, Bayesian Regularisation requires no validation
                net.divideParam.trainRatio = 85/100;
                net.divideParam.valRatio = 0/100;
                net.divideParam.testRatio = 15/100;

                % Min gradient
                % net.trainParam.min_grad = 1e-20;
                
                tic
                % GPU training does not work for Jacobian type training methods
                % st = std(y(i, :));
                % mn = mean(y(i, :));
                % trainInput = (y(i,:)-mn)/st;
                % [net, tr] = train(net, freqSample,trainInput);
                [net, tr] = train(net, freqSample, y(i, :));
                disp(['------------ Trained ',num2str(i),'/',num2str(noNetworks), ' neural networks (2x modes) (PODP) ------------'])
                timeTakenPer = toc;
                disp(['Time taken to train: ', num2str(timeTakenPer), ' seconds.'])

                pred = net(freqOut);
                % fprintf('size of pred (%d, %d)', size(pred))
                % predictions(i, :) = (pred * st) + mn;
                predictions(i, :) = pred;
                
                                
                saveLocation = sprintf('%sLossCurve_%s_l%d_n%d_m%d_Ns%d_network%d', saveFolder, solver, hiddenLayers, hiddenNeurons, noNetworks/2, size(trainSet, 2), i);
                save(saveLocation, 'tr')
                % tr.epoch
                % e = gsubtract(x,pred)
                % performance = perform(net, freqSample, pred);
                disp(['Best training MSE loss: ', num2str(tr.best_perf)])
                disp(['Best validation MSE loss: ', num2str(tr.best_vperf)])
                disp(['Best test MSE loss: ', num2str(tr.best_tperf)])
                disp(['Stopping criteria: ', tr.stop])
            end
            p = conj(cols2complex(transpose(predictions)));
            % fprintf('Size of p(w) is (%d, %d) \n', size(p))
            % new_q = H * S * p';
            new_q = H * p';
            fprintf('Size of new q(w) is (%d, %d) \n', size(new_q))
   
            timeTaken = toc(startTrain);
            disp(['Time taken to compute neural networks: ', num2str(timeTaken), ' seconds.'])
        end
    end
end