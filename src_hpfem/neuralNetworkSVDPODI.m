function [sol] = neuralNetworkSVDPODI(freqSample, H, S, G, freqOut, sizePrint, Options)
    % if loading in H, S, G
    % format long e
    % addpath(genpath('./..'))
    % load('SVDResult.mat','H','S','G');
    P = py.sys.path;
    if count(P,'src_hpfem') == 0
        insert(P,int32(0),'src_hpfem');
    end
    

    if Options.feedForwardNet==1
        y = transpose(complex2cols(conj(G)));
        noNetworks = size(y, 1);
        solver = Options.ffnSolver;
        hiddenLayers = Options.ffnLayers;
        hiddenNeurons = Options.ffnNeurons;
        predictions = zeros(noNetworks, size(freqOut, 2));
        hiddenSizes = ones(1, hiddenLayers)*hiddenNeurons
        % size(y)
        % size(freqSample)
        if Options.ffnPerMode==1
            saveFolder = sprintf('NNdata/PODI/loss_curves/loss_curves_%s_l%d_n%d_m%d_Ns%d_reg%d/', solver, hiddenLayers, hiddenNeurons, noNetworks/2, size(G, 1), 10*Options.ffnReg);
            if ~exist(saveFolder, 'dir')
                mkdir(saveFolder)
            end
            
            if Options.ffnLoggedFreq == 1
                x = log10(freqSample);
                freqOut = log10(freqOut);
            else
                x = freqSample;
            end
            
            startTrain = tic;
            for i = 1:noNetworks
                if Options.functionFitNet == 1
                    net = fitnet(hiddenSizes, solver);
                else
                    net = feedforwardnet(hiddenSizes, solver);
                end
                
                % Setup Division of Data for Training, Validation, Testing, Bayesian Regularisation requires no validation
                net.divideParam.trainRatio = Options.ffnTrainRatio;
                net.divideParam.valRatio = Options.ffnValRatio;
                net.divideParam.testRatio = Options.ffnTestRatio;

                % percentage regularisation
                net.performParam.regularization = Options.ffnReg;

                % minimum gradient
                net.trainParam.min_grad = Options.ffnMinGrad;

                % maximum epochs
                net.trainParam.epochs = Options.ffnEpochs;
                
                tic
                % GPU training does not work for Jacobian type training methods
                [net, tr] = train(net, x, y(i, :));
                disp(['------------ Trained ',num2str(i),'/',num2str(noNetworks), ' neural networks (2x modes) (PODI)------------'])
                timeTakenPer = toc;
                disp(['Time taken to train: ', num2str(timeTakenPer), ' seconds.'])

                pred = net(freqOut);
                predictions(i, :) = pred;
                
                                
                saveLocation = sprintf('%sLossCurve_%s_l%d_n%d_m%d_Ns%d_reg%d_network%d', saveFolder, solver, hiddenLayers, hiddenNeurons, noNetworks/2, size(G, 1),10*Options.ffnReg, i);
                save(saveLocation, 'tr')
                
                % e = gsubtract(x,pred)
                % performance = perform(net, freqSample, pred);
                disp(['Best training MSE loss: ', num2str(tr.best_perf)])
                disp(['Best validation MSE loss: ', num2str(tr.best_vperf)])
                disp(['Best test MSE loss: ', num2str(tr.best_tperf)])
                disp(['Stopping criteria: ', tr.stop])
            end
            new_G = conj(cols2complex(transpose(predictions)));
            saveLocationG = sprintf('%sGCurve_%s_l%d_n%d_m%d_Ns%d_reg%d', saveFolder, solver, hiddenLayers, hiddenNeurons, noNetworks/2, size(G, 1),10*Options.ffnReg);
            save(saveLocationG, 'new_G')
            sol = H * S * new_G';
            timeTaken = toc(startTrain);
            disp(['Time taken to compute neural networks: ', num2str(timeTaken), ' seconds.'])
        else
            saveFolder = sprintf('NNdata/loss_curves/PODI/loss_curves_%s_l%d_n%d_m%d_Ns%d/', solver, hiddenLayers, hiddenNeurons, noNetworks/2, size(G, 1));
            if ~exist(saveFolder, 'dir')
                mkdir(saveFolder)
            end
            tic
            net = feedforwardnet(hiddenSizes, solver);
            % Setup Division of Data for Training, Validation, Testing, Bayesian Regularisation requires no validation
            net.divideParam.trainRatio = 85/100;
            net.divideParam.valRatio = 0/100;
            net.divideParam.testRatio = 15/100;
            [net, tr] = train(net, freqSample, y);
            disp(['------------ Trained neural network for all modes------------'])
            predictions = net(freqOut);
            
            saveLocation = sprintf('%sLossCurve_%s_l%d_n%d_m%d_Ns%d_allModes', saveFolder, solver, hiddenLayers, hiddenNeurons, noNetworks/2, size(G, 1));
            save(saveLocation, 'tr')

            disp(['Best training MSE loss: ', num2str(tr.best_perf)])
            disp(['Best validation MSE loss: ', num2str(tr.best_vperf)])
            disp(['Best test MSE loss: ', num2str(tr.best_tperf)])
            disp(['Stopping criteria: ', tr.stop])


            new_G = conj(cols2complex(transpose(predictions)));
            sol = H * S * new_G';
            timeTaken = toc;
            disp(['Time taken to compute neural network: ', num2str(timeTaken), ' seconds.'])

        end
    else
        X = py.numpy.array(freqSample);

    
        G = complex2cols(conj(G));
        % y = py.numpy.array(S * G')
        
        y = py.numpy.array(G);
        % snapshots as inputs, G being the correspsonding output 

        % make log graph possibly?
        %    freqOut = log(freqOut);
        prediction_data = py.numpy.array(freqOut);
  
        % ----------------------------------------------------------------
        saveTestTrainSplit = false;
        neurons = Options.networkNeurons;
        layers = Options.networkLayers;
        solver = Options.networkSolver;
        activation = Options.networkActivation;
        segmentNo = Options.segmentNo;
        % ----------------------------------------------------------------
        neurons = cast(neurons, 'int32');
        layers = cast(layers, 'int32');
        if Options.networkSegmented==1
            segmentNo = cast(segmentNo, 'int32');
            pydata = py.funcs.model_predict_segments(X, y, prediction_data, solver, activation, neurons, layers, segmentNo);
            prediction = double(pydata);
        elseif Options.networkPerMode==1
            pydata = py.funcs.model_predict_per_mode(X, y, prediction_data, solver, activation, neurons, layers);
            prediction = double(pydata);
        else
            pydata = py.funcs.model_predict(X, y, prediction_data, solver, activation, neurons, layers, saveTestTrainSplit);
            prediction = double(pydata{1})
        end     
        %    fprintf("size of prediction is: (%d, %d) \n", size(prediction))
        new_G = conj(cols2complex(prediction));
        % X_train = double(pydata{3});
        sol = H * S * new_G';
        if sizePrint == true
            fprintf("size of H: (%d, %d) \n", size(H))
            fprintf("size of S: (%d, %d) \n", size(S))
            fprintf("size of G: (%d, %d) \n", size(G))
            fprintf("size of snapshots: (%d, %d) \n", size(freqSample))
            fprintf("size of G as [real, complex]: (%d, %d) \n", size(G))
            fprintf("size of X_train is: (%d, %d) \n", size(X_train))
            fprintf("size of new G hermitian is: (%d, %d) \n", size(new_G'))
            fprintf("size of sol is: (%d, %d) \n", size(sol))
        end
    end
end
                             