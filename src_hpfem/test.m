load('../data/SVDResult_Ns180_m20.mat', 'H', 'S', 'G', 'trainPODPNN', 'trainPODPNNFULL')



% Neural network
neuralNetworkPODI    = 0;      % Consider a neural network for the SVD in PODI interpolating function (1) or lagrangian (0)
neuralNetworkPODP    = 1;      % Consider a neural network for the SVD in PODP (1) or normal projection (0)

feedForwardNet   = 1;          % Consider feedforwardnet for the neural network implementation (1) or scikit-learn regressor (0)
functionFitNet   = 0;          % Consider a function fitting neural network implementation (1) or ffn (0)
ffnLayers        = 2;          % Number of hidden layers
ffnNeurons       = 16;         % Number of hidden neurons
ffnSolver        = 'trainbr';  % Training function used for back propagation
ffnPerMode       = 1;          % Train a NN per mode (1) or just one (0

Options.neuralNetworkPODI=neuralNetworkPODI;
Options.neuralNetworkPODP=neuralNetworkPODP;
Options.feedForwardNet=feedForwardNet;
Options.functionFitNet=functionFitNet;
Options.ffnNeurons=ffnNeurons;
Options.ffnSolver=ffnSolver;
Options.ffnLayers=ffnLayers;
Options.ffnPerMode=ffnPerMode;



del_f1 = 10;
del_f2 = 50;
s1 = 10:del_f1:1000;
s2 = 1000:del_f2:5000;
freqSample = [s1, s2(1:end-1)];
fprintf("Number of snapshots: (%d, %d) \n", size(freqSample))

del_fout = 10;
N_o = 499;
freqOut = linspace(15,15+(N_o-1)*del_fout, N_o);
disp("orthogonality check")
H(:, 1) * H(:, 2)' - eye(size(H,2))
disp("(H^m)^* (H^m) - I")
H'*H - eye(size(H,2))



% A = transpose(complex2cols(transpose(trainPODPNN)));
% diff1 = trainPODPNN - (S * G');
% disp("(H^m)'D - (S^m)(G^m)*")
% mean(diff1, 'all')
% disp("(H)'D - (S^m)(G^m)*")
% diff2 = trainPODPNNFULL - (S * G')

% A(11, :)
% B = transpose(cols2complex(transpose(A)));
% B == trainPODPNN;



% neuralNetworkSVDPODP(freqSample, H, S, G,trainPODPNN, freqOut, Options);





