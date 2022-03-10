addpath(genpath('./..'))
% N_s value to load and interpolate SVD
N_s = 1296
m = 20

load(['data/SVDResult_Ns',num2str(N_s),'_m',num2str(m)],'H', 'S', 'G');
oldG = G;
if N_s == 1296
    del_f1 = 2;
    del_f2 = 5;
elseif N_s == 599
    del_f1 = 5;
    del_f2 = 10;
else
    del_f1 = 1;
    del_f2 = 1;
end

s1 = 10:del_f1:1000;
s2 = 1000:del_f2:5000;
snapshots = [s1, s2(1:end-1)];

% del_f1_o = 1;
% del_f2_o = 3;

% s1_o = 10:del_f1_o:1000;
% s2_o = 1000:del_f2_o:5000;
% freqout = [s1_o, s2_o(1:end-1)];
% N_s_new = length(freqout);

% freqout = double.empty

% for i=1:length(snapshots)-1
%     freqout(end+1) = ceil((snapshots(i+1) + snapshots(i))/2);
% end
% freqout
N_o = 2500;
freqout = linspace(10, 5000, N_o);
N_s_new = length(freqout);


G = LagrangeInterpolationG(snapshots, H, S, oldG, freqout)
fprintf("Number of old snapshots: %d \n", N_s)
fprintf("size of old G is: (%d, %d) \n", size(oldG))
fprintf("Number of new snapshots: %d \n", N_s_new)
fprintf("size of new G is: (%d, %d) \n", size(G))


save(['../data/SVDResult_Ns', num2str(N_s_new),'_m',num2str(m),'_int'],'H','S','G');



