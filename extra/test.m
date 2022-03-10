

% N_s = 180
% del_f1 = 5;
% del_f2 = 10;
% calc = (1000-10)/del_f1 + (5000-1000)/del_f2
% while calc > N_s
%     del_f1 = del_f1 + 1;
%     del_f2 = del_f2 + 1;
%     calc = (1000-10)/del_f1 + (5000-1000)/del_f2
% end
% disp(del_f1)
% disp(del_f2)


% ns 2324 1296 599 359 180 90 45 23
% delf1 1 2 5 5 10 20 40 80
% delf2 3 5 10 25 50 100 200 400

del_f1 = 10;
del_f2 = 50;
s1 = 10:del_f1:1000;
s2 = 1000:del_f2:5000;
sample = [s1, s2(1:end-1)];

sample
length(sample)

% l = 5;
% n = 5;

% s = ones(1, l)*n

% tol = 1e-14
% disp(['tol_',num2str(tol)])

% freqSample = 5:5:5432;
% for j=1:length(freqSample)
%     freq = freqSample(j);
%     if rem(j, round(length(freqSample)/10)) == 0
%         disp(['---------------- Linear system solved for ', num2str(j),'/',num2str(length(freqSample)),' samples ----------------'])
%     end
% end


Z = [0.5i 43; 1+3i 4i; -2.2 3i];
%[complex2cols(Z(:,1)) complex2cols(Z(:,2))]
% complex2cols(Z)
% cols2complex(complex2cols(Z))
 