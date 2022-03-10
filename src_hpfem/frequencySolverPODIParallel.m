function [Dynamic]=frequencySolverPODIParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,Damp,dampRatio,CondFactorOut,CondFactorChoice,CondFactorSample,freqSample,nModes,Ncores)
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
    DynamicSampleTotal=cell(nCond);

    % Start parallel pool
    parpool(Ncores)
    % Copy matrices to workers for parallel efficiency
    MP=parallel.pool.Constant(M);
    KP=parallel.pool.Constant(K);
    ResidP=parallel.pool.Constant(Resid);
    CRegP=parallel.pool.Constant(CReg);
    CpreRegP=parallel.pool.Constant(CpreReg);
    clear CReg CpreReg

    for i=1:nCond
    DynamicSample = zeros(lenEM+ndirEM,length(freqSample));
    C=CondFactorSample(i,1)*Ccond{1};
    Cpre=CondFactorSample(i,1)*Ccondpre{1};
    for nBody=2:NmechBodies
    C=C+CondFactorSample(i,nBody)*Ccond{nBody};
    Cpre=Cpre+CondFactorSample(i,nBody)*Ccondpre{nBody};
    end
    CP=parallel.pool.Constant(C);
    CpreP=parallel.pool.Constant(Cpre);
    q=parallel.pool.DataQueue;
    afterEach(q, @disp)

        parfor j=1:length(freqSample)
            
            % Define angular frequency (rad/s) from frequency (Hz)
            freq=freqSample(j);
            omega       = freq*2*pi;
            % Compute the K, C and M weights
            w0         = 1;
            w1         = complex(0,1)*omega;
            w2         = -omega^2;
            

            
            [U] = linearSystemSolverEM(Unknown,ProblemData,w0,w1,w2,...
                MP.Value,CP.Value,CRegP.Value,KP.Value,ResidP.Value,CpreP.Value,CpreRegP.Value,Options);
            disp(['---------------- Linear system solved for ',num2str(j),'/',num2str(length(freqSample)),' samples ----------------'])


            if Options.Non0Dir==1
                % Include Dirichlet Values in the Solution
                probFlag=1;
                [AAC] = initialGuessP(Mesh,Basis,Quadrature,Unknown.EM,ProblemData,omega,probFlag);
                
                %-------------------------------------------------------------------------
                % Initial guess of the mechanical field
                %--------------------------------------------------------------------------
                
                % Define the problem flag
                probFlag=2;
                
                % Initial guess for the mechanical problem
                [UAC] = initialGuessP(Mesh,Basis,Quadrature,Unknown.Mech,ProblemData,omega,probFlag);
                X2=[AAC;UAC];
                U=U+X2(1:lenEM+ndirEM);
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
            DynamicSample(:,j)= U;
            
        end
        DynamicSampleTotal{i}=DynamicSample;
    end
    clear CP CRegP CpreRegP CpreP
    % Rearrange solution into matrix
    DynamicSample=zeros(lenEM+ndirEM,nCond*length(freqSample));
    for i=1:nCond
        DynamicSample(:,(i-1)*length(freqSample)+1:i*length(freqSample))=DynamicSampleTotal{i};
    end
    %====================================================================================
    % Compute TSVD to obtain POD modes
    %====================================================================================
    
    % default tolerence: 1e-14
    tol = 1e-14;
    [H,S,G]=svds(DynamicSample(1:lenEM,:),nModes, 'largest', 'Tolerance', tol);
    %    trainPODPNN = H'*DynamicSample;

    % Save the result of the SVD
    if Options.svdSave == 1
        save(['data/',Options.svdSaveName],'H','S','G');
    else        
        if tol == 1e-14
            if Options.chooseSample == 'marcos'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes)],'H','S','G','trainPODPNN');
            elseif Options.chooseSample == 'log spaced'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz_logspace'],'H','S','G');
            end
        else
            if Options.chooseSample == 'marcos'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes), '_tol',num2str(tol)],'H','S','G');
            elseif Options.chooseSample == 'log spaced'
                save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes),'_tol',num2str(tol),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz_logspace'],'H','S','G');
            end
        end
    end
    clear DynamicSample
end
clear C Cpre Ccondpre 
disp('Offline stage complete.')
%====================================================================================
% Online Stage
%====================================================================================
if Options.offlineOnly == 1
    Dynamic = double.empty;
    disp('No online stage to be computed. Program terminated.')
else
    disp('------------- Online stage ------------')
    HP=parallel.pool.Constant(H);
    SP=parallel.pool.Constant(S);
    GP=parallel.pool.Constant(G);
    % clear H
    [nCond,~]=size(CondFactorOut);
    DynamicTotal=cell(nCond);

    stiffRed=H'*K(1:lenEM,1:lenEM)*H;
    ResidRed=H'*Resid(1:lenEM);
    
    stiffRedP=parallel.pool.Constant(stiffRed);
    ResidRedP=parallel.pool.Constant(ResidRed);
    clear stiffRed
    clear ResidRed


    NmechBodies=ProblemData.NmechBodies;
    for i=1:nCond
        Dynamic = zeros(nTotal,length(freqOut));
        parfor j=1:length(freqOut)
            
            % Define angular frequency (rad/s) from frequency (Hz)
            freq=freqOut(j);
            omega       = freq*2*pi;
            % Compute the K, C and M weights
            w0         = 1;
            w1         = complex(0,1)*omega;
            w2         = -omega^2;

            AredWhole=LagrangeInterpolation(freqSample,HP.Value,SP.Value,GP.Value,freqOut);
            %          AFull=linearSystemSolverEM(Unknown,ProblemData,w0,w1,CP.Value,CRegP.Value,KP.Value,ResidP.Value,CpreP.Value,CpreRegP.Value);
            % %         %DynamicEMFull(:,j)=AFull;
            %          ErrorVector(j)=norm(AredWhole-AFull(1:lenEM))/norm(AFull(1:lenEM));
            
            if Options.SplitMech==1
                [U] = linearSystemSolverMech(Unknown,AredWhole(:,j),w0,w1,w2,...
                                             MP.Value,KP.Value,ResidP.Value,dampRatio,NmechBodies);
            else
                [U] = linearSystemSolverMechNS(Unknown,AredWhole(:,j),w0,w1,w2,...
                                               MP.Value,KP.Value,ResidP.Value,dampRatio);  
            end
            if Options.Non0Dir==1
                % Include Dirichlet Values in the Solution
                probFlag=1;
                [AAC] = initialGuessP(Mesh,Basis,Quadrature,Unknown.EM,ProblemData,omega,probFlag);
                
                %-------------------------------------------------------------------------
                % Initial guess of the mechanical field
                %--------------------------------------------------------------------------
                
                % Define the problem flag
                probFlag=2;
                
                % Initial guess for the mechanical problem
                [UAC] = initialGuessP(Mesh,Basis,Quadrature,Unknown.Mech,ProblemData,omega,probFlag);
                X2=[AAC;UAC];
                U=U+X2;
            end 
            %---------------------------------------------------------------------
            % Store the solution at each frequency
            %---------------------------------------------------------------------
            % Adjust the Frequency domain solutions
            Dynamic(:,j)= U;        
        end
        DynamicTotal{i}=Dynamic;
    end
    % Rearrange solution into matrix
    Dynamic=zeros(nTotal,nCond*length(freqOut));
    for i=1:nCond
        Dynamic(:,(i-1)*length(freqOut)+1:i*length(freqOut))=DynamicTotal{i};
    end
end
clear H U Ared AredWhole dampRed stiffRed ResidRed
clear M K Resid C COVC C77K C4K Cpre CpreReg CReg Cpre77K Cpre4K CpreOVC

