clear all
clc

files = dir('*.mat');

plotPower = 1;      % Plot power plots (1) or not (0)
% plotCurl = 0;       % Plot norm of curl A plots (1) or not (0)

maxF = 5000;        % Max frequency

% Initialise import indexing
powerCount = 1;
% curlCount = 1;

% Import data
for i = 1:length(files)
    load(files(i).name);
    if contains(convertCharsToStrings(files(i).name),'power')
%         if plotPower == 1
            power(powerCount) = IntegratedFields;
            powerCount = powerCount + 1;
%         end
%     elseif contains(convertCharsToStrings(files(i).name),'norm')
%         if plotCurl == 1
%             curl(curlCount) = IntegratedNormCurlA;
%             curlCount = curlCount + 1;
%         end
    end
end

% Import power data if flag selected

if plotPower == 1
    % Initialise empty arrays
    outPower4K = zeros(length(power(1).OutPower4K),length(power));
    outPower77K = zeros(length(power(1).OutPower77K),length(power));
    outPowerOVC = zeros(length(power(1).OutPowerOVC),length(power));
    dispNorm4K = zeros(length(power(1).DisplacementNorm4K),length(power));
    dispNorm77K = zeros(length(power(1).DisplacementNorm77K),length(power));
    dispNormOVC = zeros(length(power(1).DisplacementNormOVC),length(power));
    
    % Populate empty arrays
    for i = 1:length(power)
        outPower4K(:,i) = power(i).OutPower4K;
        outPower77K(:,i) = power(i).OutPower77K;
        outPowerOVC(:,i) = power(i).OutPowerOVC;
        dispNorm4K(:,i) = power(i).DisplacementNorm4K;
        dispNorm77K(:,i) = power(i).DisplacementNorm77K;
        dispNormOVC(:,i) = (power(i).DisplacementNormOVC);
    end
end

% % Import curl data if flag selected
% 
% if plotCurl == 1
%     % Initialise empty arrays
%     outNormCurlA4K = zeros(length(curl(1).OutNormCurlA4K),length(curl));
%     outNormCurlA44K = zeros(length(curl(1).OutNormCurlA77K),length(curl));
%     outNormCurlAOVC = zeros(length(curl(1).OutNormCurlAOVC),length(curl));
%     
%     % Populate empty arrays
%     for i = 1:length(curl)
%         outNormCurlA4K(:,i) = curl(i).OutNormCurlA4K;
%         outNormCurlA44K(:,i) = curl(i).OutNormCurlA77K;
%         outNormCurlAOVC(:,i) = curl(i).OutNormCurlAOVC;
%     end
% end

% End of data import

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kinetic Energy calculation

KE4K = dispNorm4K.^2;
KE77K = dispNorm77K.^2;
KEOVC = dispNormOVC.^2;

df = maxF / (length(dispNorm4K) + 1);
f = 0;

for i = 1:length(outPower4K)
    f = f + df;
    KE4K(i,:) = 0.5 * 4 * 2710 * ((2 * pi() * f)^2) * (KE4K(i,:));
    KE77K(i,:) = 0.5 * 4 * 22698 * ((2 * pi() * f)^2) * (KE77K(i,:));
    KEOVC(i,:) = 0.5 * 4 * 27900 * ((2 * pi() * f)^2) * (KEOVC(i,:));
end


% Results plotting

xPower = linspace(0,maxF,length(outPower4K));
% xCurl = linspace(0,5000,length(outPower4K));

figure

% dispLegend = {'disp = 0.100m', 'disp = 0.080m', 'disp = 0.020m', 'disp = 0.005m', 'disp = 0.001m', 'disp = 0.000m'};
% dispLegend = {'disp = 0.000m', 'disp = 0.0005m', 'disp = 0.001m', 'disp = 0.005m', 'disp = 0.020m', 'disp = 0.080m', 'disp = 0.100m'};
dispLegend = {'disp = 0.0000m', 'disp = 0.0005m', 'disp = 0.0010m', 'disp = 0.0050m', 'disp = 0.0200m'};

t1 = tiledlayout(3,2);

ax1 = nexttile(5);
semilogy(ax1,xPower,outPower4K,'LineWidth',2.0)
title('Dissipated Power 4K')
xlabel('$Frequency \; (Hz)$','Interpreter','Latex')
ylabel('$P^0 (W)$','Interpreter','Latex')
ylim([10e-3 10e9])
grid on
legend(dispLegend{:}, 'location', 'nw');


ax2 = nexttile(3);
semilogy(ax2,xPower,outPower77K,'LineWidth',2.0)
title('Dissipated Power 77K')
xlabel('$Frequency \; (Hz)$','Interpreter','Latex')
ylabel('$P^0 (W)$','Interpreter','Latex')
ylim([10e-1 10e10])
grid on
legend(dispLegend{:}, 'location', 'nw');

ax3 = nexttile(1);
semilogy(ax3,xPower,outPowerOVC,'LineWidth',2.0)
title('Dissipated Power OVC')
xlabel('$Frequency \; (Hz)$','Interpreter','Latex')
ylabel('$P^0 (W)$','Interpreter','Latex')
ylim([10e-2 10e7])
grid on
legend(dispLegend{:}, 'location', 'nw');

ax4 = nexttile(6);
semilogy(ax4,xPower,KE4K,'LineWidth',2.0)
title('Kinetic Energy 4K')
xlabel('$Frequency (Hz)$','Interpreter','Latex')
ylabel('$Kinetic \; Energy \; (J)$','Interpreter','Latex')
ylim([10e-10 10e14])
grid on
legend(dispLegend{:}, 'location', 'nw');

ax5 = nexttile(4);
semilogy(ax5,xPower,KE77K,'LineWidth',2.0)
title('Kinetic Energy 77K')
xlabel('$Frequency (Hz)$','Interpreter','Latex')
ylabel('$Kinetic \; Energy \; (J)$','Interpreter','Latex')
ylim([10e-9 10e15])
grid on
legend(dispLegend{:}, 'location', 'nw');

ax6 = nexttile(2);
semilogy(ax6,xPower,KEOVC,'LineWidth',2.0)
title('Kinetic Energy OVC')
xlabel('$Frequency (Hz)$','Interpreter','Latex')
ylabel('$Kinetic \; Energy \; (J)$','Interpreter','Latex')
ylim([10e-12 10e16])
set(gcf, 'Position',  [0, 0, 800, 800])
grid on
legend(dispLegend{:}, 'location', 'nw');