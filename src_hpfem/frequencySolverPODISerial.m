function [Dynamic]=frequencySolverPODISerial(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,Damp,dampRatio,CondFactorOut,CondFactorChoice,CondFactorSample,freqSample,nModes)
% Function used to initialise the linearised solver of the coupled equation
% set

disp('--------------------------------------------')
disp('Running the Frequency Domain solver')
disp('--------------------------------------------')


%=========================================================================
% Import the data from the data structures
%=========================================================================
couple        = Options.Couple;
Assemble      = Options.Assemble;
SourceMapping = Options.SourceMapping;
freqSweep     = Options.freqSweep;
NmechBodies   = ProblemData.NmechBodies;

% Initialize the unknown matrices
M      = [];
CReg   = [];
K      = [];


%=========================================================================
% Solve the transient fields
%=========================================================================
% Initialise the transient analysis flag
ProblemData.probstatic = 0;
% Store angular frequency in ProblemData structure (used for BCs when solving for 1 frequency only)
ProblemData.matr.omega=2*pi*freqOut(1);
% Compute the system matrices for the first frequency
lenEM=Unknown.EM.nunkt;
lenM=Unknown.Mech.nunkt;
ndirEM=Unknown.EM.npec;
ndirMech=Unknown.Mech.npec;
nTotal=ndirEM+ndirMech+lenEM+lenM;


% Store exact spase matrix vectors size in Unknown structure
Unknown.system.nSparse=Static.nSparse;

%-------------------------------------------------------------------------
% Frequency sweep
%-------------------------------------------------------------------------

%---------------------------------------------------------------------
% Linear system solver
%---------------------------------------------------------------------
% Include Dirichlet values in solution
if Options.Non0Dir==1 && freqSweep==0
    % Initial guess of the solution field
    ProblemData.ProbFlag=1;
    [AAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.EM,ProblemData,2*pi*freqOut);
    % Define the problem flag
    ProblemData.probFlag=2;
    
    % Initial guess for the mechanical problem
    [UAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.Mech,ProblemData,2*pi*freqOut);
    X=[AAC;UAC];
else
    X=zeros(nTotal,1);
end

%=========================================================================
% System assembly
%=========================================================================
fileName=sprintf('MatricesDynamic_%d',ProblemData.order);
fileName=strcat(fileName,ProblemData.jb.job);
if Assemble==1
    disp('Constructing the Stiffness and Forcing Matrices')
    %[M1,C1,Cpre1,CReg1,CpreReg1,K1,Resid1,~,~,~,~,~]=elementLoop4(Static,Mesh,Basis,Quadrature,Unknown,UnknownStatic,ProblemData,X,couple,freqSweep);
    % Assemble the system components into the tangent matrix and Residual
    [M,Ccond,Ccondpre,CReg,CpreReg,K,Resid,~,~,~,~,~]=elementLoop(Static,StaticCurrent,Mesh,Basis,Quadrature,Unknown,UnknownCurrent,UnknownStatic,ProblemData,X,couple,freqSweep,SourceMapping);
    save(['data/',fileName],'M','Ccond','Ccondpre','CReg','CpreReg','K','Resid');
else
    load(['data/',fileName],'M','Ccond','Ccondpre','CReg','CpreReg','K','Resid');
    %Resid=Res;
end

if Options.Offline==0
    if Options.svdSave==1
        load(['data/', Options.svdSaveName], 'H', 'S', 'G');
    else
        if Options.chooseSample == 'marcos'
            load(['data/SVDResult_Ns',num2str(size(freqSample, 2)),'_m',num2str(nModes)],'H', 'S', 'G');
        elseif Options.chooseSample == 'log spaced'
            load(['data/SVDResult_Ns', num2str(size(freqSample, 2)),'_m',num2str(nModes),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz_logspace'], 'H', 'S', 'G');
        end
    end
else
    [nCond,~]=size(CondFactorSample);
    DynamicSample = zeros(lenEM+ndirEM,nCond*length(freqSample));
    count=0;
    for i=1:nCond
    C=CondFactorOut(i,1)*Ccond{1};
    Cpre=CondFactorOut(i,1)*Ccondpre{1};
    for nBody=2:NmechBodies
    C=C+CondFactorOut(i,nBody)*Ccond{nBody};
    Cpre=Cpre+CondFactorOut(i,nBody)*Ccondpre{nBody};
    end
        for j=1:length(freqSample)
            
            count=count+1;
            % Define angular frequency (rad/s) from frequency (Hz)
            freq=freqSample(j);
            omega       = freq*2*pi;
            % Compute the K, C and M weights
            w0         = 1;
            w1         = complex(0,1)*omega;
            w2         = -omega^2;
            
            if Options.Non0Dir==1
                % Include Dirichlet Values in the Solution
                ProblemData.ProbFlag=1;
                [AAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.EM,ProblemData,omega);
                
                %-------------------------------------------------------------------------
                % Initial guess of the mechanical field
                %--------------------------------------------------------------------------
                
                % Define the problem flag
                ProblemData.probFlag=2;
                
                % Initial guess for the mechanical problem
                [UAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.Mech,ProblemData,omega);
                X=[AAC;UAC];
            end
            
            [U] = linearSystemSolverEM(Unknown,ProblemData,w0,w1,w2,...
                M,C,CReg,K,Resid,Cpre,CpreReg,Options);

            U=U+X(1:lenEM+ndirEM);
            if rem(j, round(length(freqSample)/10)) == 0
                disp(['---------------- Linear system solved for ', num2str(j),'/',num2str(length(freqSample)),' samples ----------------'])
            end

            
            %---------------------------------------------------------------------
            % Post Processing - compute output power and kinetic energy of
            %                   conducting components
            %---------------------------------------------------------------------
            
            % Compute the integrated values in time
            %     if fieldCalc == 1
            %         Dynamic{j}.Field =integratedFieldOutput(Matrices,Unknown,Options,X,dXdt,dXdt2);
            %     end
            %
            
            %---------------------------------------------------------------------
            % Store the solution at each frequency
            %---------------------------------------------------------------------
            % Adjust the Frequency domain solutions
            DynamicSample(:,count)= U;
            
        end
    end
    %====================================================================================
    % Compute TSVD to obtain POD modes
    %====================================================================================
        
    % default tolerence: 1e-14
    tol = 1e-14;
    [H,S,G]=svds(DynamicSample(1:lenEM,:),nModes, 'largest', 'Tolerance', tol);
    trainPODPNN = H'*DynamicSample;
    % Save the result of the SVD
        % Save the result of the SVD
    if Options.svdSave == 1
        save(['data/',Options.svdSaveName],'H','S','G');
    else        
        if tol == 1e-14
            if Options.chosenSample == 'marcos'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes)],'H','S','G','trainPODPNN');
            elseif Options.chosenSample == 'log spaced'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz_logspace'],'H','S','G');
            end
        else
            if Options.chosenSample == 'marcos'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes), '_tol',num2str(tol)],'H','S','G');
            elseif Options.chosenSample == 'log spaced'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes),'_tol',num2str(tol),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz_logspace'],'H','S','G');

            end
        end
    end

clear DynamicSample
end

clear Cpre Ccondpre 

%====================================================================================
% Online Stage
%====================================================================================
H=H(1:lenEM,:);
[nCond,nBodies]=size(CondFactorOut);


Dynamic = zeros(nTotal,nCond*length(freqOut));

if Options.offlineOnly==1
    disp('No online stage to be computed. Program terminated.')    
else
    
    NmechBodies=ProblemData.NmechBodies;
    count=0;
    tStartLoops=tic;

    
    % neural network solution or just approxD
    if Options.approxD==1
        AredWhole=H*S*G';
    elseif Options.neuralNetworkPODI==1
        prints = false;
        AredWhole=neuralNetworkSVDPODI(freqSample,H,S,G,freqOut, prints, Options);
    else
        AredWhole=LagrangeInterpolation(freqSample,H,S,G,freqOut, Options);
        % AredWhole=LinearInterpolation(freqSample,H,S,G,freqOut);
    end
    if Options.networkOnly==1
        disp('Only neural networks trained with metrics saved. Program terminated.')
    else
        for i=1:nCond
            
            for j=1:length(freqOut)
                
                count=count+1;
                % Define angular frequency (rad/s) from frequency (Hz)
                freq=freqOut(j);
                omega       = freq*2*pi;
                % Compute the K, C and M weights
                w0         = 1;
                w1         = complex(0,1)*omega;
                w2         = -omega^2;
                if Options.Non0Dir==1
                    % Include Dirichlet Values in the Solution
                    ProblemData.ProbFlag=1;
                    [AAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.EM,ProblemData,omega);
                    
                    %-------------------------------------------------------------------------
                    % Initial guess of the mechanical field
                    %--------------------------------------------------------------------------
                    
                    % Define the problem flag
                    ProblemData.probFlag=2;
                    
                    % Initial guess for the mechanical problem
                    [UAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.Mech,ProblemData,omega);
                    X=[AAC;UAC];
                end 
                
                if Options.SplitMech==1
                    [U] = linearSystemSolverMech(Unknown,AredWhole(:,j),w0,w1,w2,...
                                                 M,K,Resid,dampRatio,NmechBodies);
                else
                    [U] = linearSystemSolverMechNS(Unknown,AredWhole(:,j),w0,w1,w2,...
                                                   M,K,Resid,dampRatio);  
                end
                U=U+X;
                %---------------------------------------------------------------------
                % Store the solution at each frequency
                %---------------------------------------------------------------------
                % Adjust the Frequency domain solutions
                Dynamic(:,count)= U;
                disp(['------------ Conducting factor choice ',num2str(i),': ',num2str(j),'/',num2str(length(freqOut)), ' frequencies solved. ------------'])
            end
        end
    end
end
clear H U Ared AredWhole dampRed stiffRed ResidRed
clear M K Resid C COVC C77K C4K Cpre CpreReg CReg Cpre77K Cpre4K CpreOVC


