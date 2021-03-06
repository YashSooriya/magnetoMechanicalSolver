function postStaticSolver(dirDispNew, dampRatio, orderEM, orderMech,method)

delete(gcp('nocreate'))

load('temp/postStaticSolverData.mat');

ProblemData.non0=dirDispNew;
Options.dirDisp=dirDispNew;


ticInit = tic;
toccount = 1;
Toc = [];

% Dynamic Solver (Frequency domain)
if POD==1
    if PODP==1
        [Dynamic]=frequencySolverPODPParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,dampChoice,dampRatio,CondFactorOut,CondFactorChoice,CondFactorSample,freqSample,nModes,Ncores);
    else
        [Dynamic]=frequencySolverPODIParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,dampChoice,dampRatio,CondFactorOut,CondFactorChoice,CondFactorSample,freqSample,nModes,Ncores);
    end
else
    [Dynamic, Toc, toccount]=frequencySolverFullParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,dampChoice,dampRatio,CondFactorOut,CondFactorChoice,Ncores, ticInit, Toc, toccount);
end

currDate = replace(datestr(datetime), {':',' '}, {'_','.'});
folder = ['data/dynamicData/',currDate,'/'];
mkdir(folder)
writetable(struct2table(Options),[folder,'Options.txt'])
saveFile=[folder,'postDynamicSolverData'];
save(saveFile);
disp(['Saved to ', saveFile])

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%==========================================================================
% Post processing
%==========================================================================
clear StaticCurrent UnknownCurrent
% Extract static solution
solStatic=Static.sol(:,1);

if  linePlotOn==1
    % Line Plots
    myplot3(Mesh,Unknown,ProblemData,Dynamic(:,1),problem,0,freqRange)
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%=========================================================================
% ParaView (3D plots)
%=========================================================================
if paraview==1
    % 3D plots for the dynamic solution (Paraview)
    if FixedPoint==1
        pointvalue_multi(Mesh,Unknown,ProblemData,Dynamic(:,1),0,freqRange);
    else
        pointvalue_multi_Staggered(Mesh,Unknown,ProblemData,Dynamic(:,1),0,freqOut,UnknownStatic,solStatic);
      %  pointvalue_multi_Staggered_Air(Mesh,Unknown,ProblemData,Dynamic(:,1),0,freqOut);
    end
    % 3D plots for the static solution (Paraview)
    pointvalue_multi(Mesh,UnknownStatic,ProblemData,solStatic,1,freqOut);
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%==========================================================================




%=========================================================================
% Compute error of solution (Full order) if analytical solution is available
%=========================================================================
if errorsOn==1
     hcurlerr=zeros(length(freqOut),1);
     [L2_error_mech,Energy_error_mech,SNS_error_mech] =ErrorStatic(Mesh,UnknownStatic,Basis,Quadrature,solStatic,ProblemData,1,50);
    for j=1:length(freqOut)
        freq=freqOut(j);
        sol=Dynamic(:,j);
        [hcurlerr(j),DisplacementNorm(j)] =ErrorDynamic(Mesh,Unknown,Basis,Quadrature,sol,ProblemData,0,freq);
    end
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%=========================================================================

%=========================================================================
% Compute dissipated power and kinetic Energy
%=========================================================================
    if fieldCalc==1
    disp('Calcuating Power and Energy ...')
    [IntegratedFields]= PowerAndEnergyComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic,method);
    currDate = replace(datestr(datetime), {':',' '}, {'_','.'});
    folder = ['data/powerEnergy/',currDate,'/'];
    mkdir(folder)
    writetable(struct2table(Options),[folder,'Options.txt'])
    saveFile=[folder,'FrequencySweepMHIGradXPowerEnergy'];
    save(saveFile,'IntegratedFields');
    disp(['Saved to ', saveFile])
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%=========================================================================
% Compute the Norm of A field
%=========================================================================

if normA==1
    disp('Calcuating Norm A ...')
    [IntegratedNormA]= NormAComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic);
    currDate = replace(datestr(datetime), {':',' '}, {'_','.'});
    folder = ['data/normA/',currDate,'/'];
    mkdir(folder)
    writetable(struct2table(Options),[folder,'Options.txt'])
    saveFile=[folder,'FrequencySweepMHIGradXNormA'];
    save(saveFile,'IntegratedNormA');
    disp(['Saved to ', saveFile])
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%=========================================================================

%=========================================================================
% Compute the Norm of the curlA field
%=========================================================================

if normCurlA==1
    disp('Calcuating Norm curlA ...')
    [IntegratedNormCurlA]= CurlNormAComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic);
    currDate = replace(datestr(datetime), {':',' '}, {'_','.'});
    folder = ['data/normCurlA/',currDate,'/'];
    mkdir(folder)
    writetable(struct2table(Options),[folder,'Options.txt'])
    saveFile=[folder,'FrequencySweepMHIGradXNormCurlA'];
    save(saveFile,'IntegratedNormCurlA');
    disp(['Saved to ', saveFile])
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%=========================================================================

%=========================================================================
% Compute magnetic vector potential
%=========================================================================
if normA==1
    [IntegratedNormA]= NormAComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic);
    fileName=sprintf('FrequencySweepParallelFullTestNormA%dto%dHz_damp%s_q%d_p%d_CondFactor%sOVCNS',freqOut(1),freqOut(length(freqOut)),dampChoice,orderEM,orderMech,CondFactorChoice);
    save(fileName,'IntegratedNormA');
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%=========================================================================
% Compute error of POD solution with respect to full order (e_2)
%=========================================================================
if ErrorPOD==1
DynamicPOD=Dynamic;
fileName=sprintf('FullOrderSolution');
load(fileName,'Dynamic');
[nCond,~]=size(CondFactorOut);
lenEM=Unknown.EM.nunkt;
lenM=Unknown.Mech.nunkt;
ndirEM=Unknown.EM.npec;
ndirMech=Unknown.Mech.npec;
nTotal=ndirEM+ndirMech+lenEM+lenM;
ErrorPOD=zeros(nTotal,nCond*length(freqOut));
ErrorPODEM=zeros(lenEM,nCond*length(freqOut));
ErrorPODMech=zeros(lenM,nCond*length(freqOut));
% Set counter to zero
count=0;
for i=1:nCond
    for j=1:length(freqOut)
        count=count+1;
        ErrorPOD(:,count)=norm(DynamicPOD(:,count)-Dynamic(:,count))/norm(DynamicPOD(:,count));
        ErrorPODEM(:,count)=norm(DynamicPOD(1:lenEM,count)-Dynamic(1:lenEM,count))/norm(DynamicPOD(1:lenEM,count));
        ErrorPODMech(:,count)=norm(DynamicPOD(lenEM+ndirEM+1:lenEM+ndirEM+lenM,count)-Dynamic(lenEM+ndirEM+1:lenEM+ndirEM+lenM,count))/norm(DynamicPOD(lenEM+ndirEM+1:lenEM+ndirEM+lenM,count));
    end
end
fileName=sprintf('ErrorPODCase1');
save(fileName,'ErrorPOD','ErrorPODEM','ErrorPODMech');
end

Toc(toccount) = toc(ticInit);

Toc = Toc.';
currDate = replace(datestr(datetime), {':',' '}, {'_','.'});
folder = ['data/runTime/',currDate,'/'];
mkdir(folder)
writetable(array2table(Toc),[folder,'ticktoc.csv'])