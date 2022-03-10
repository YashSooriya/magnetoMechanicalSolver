function [Dynamic]=frequencySolverPODPParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,Damp,dampRatio,CondFactorOut,CondFactorChoice,CondFactorSample,freqSample,nModes,Ncores)
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
    %save(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','Kd','Cd','CdReg','Md','-v7.3');
else
    load(['data/',fileName],'M','Ccond','Ccondpre','CReg','CpreReg','K','Resid');

    % load(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','Kd','Cd','CdReg','Md');
    %Resid=Res;
end

if Options.Offline==0
    if Options.svdSave==1
        load(['data/', Options.svdSaveName], 'H', 'S', 'G', 'trainPODPNN');
    else
        if freqSample(length(freqSample)) > 4500
            load(['data/SVDResult_Ns',num2str(size(freqSample, 2)),'_m',num2str(nModes)],'H', 'S', 'G', 'trainPODPNN');
        else
            load(['data/SVDResult_Ns', num2str(size(freqSample, 2)),'_m',num2str(nModes),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz'], 'H', 'S', 'G', 'trainPODPNN');
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
    clear CpreReg
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
    [H,S,G]=svds(DynamicSample,nModes,'largest','Tolerance', tol);
    % [tall(Hfull), tall(Sfull), tall(Gfull)] = svd(DynamicSample);
    trainPODPNN = H'*DynamicSample;
    % trainPODPNNFULL = Hfull'*DynamicSample;

    % Save the result of the SVD
    if tol == 1e-14
        save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz'],'H','S','G', 'trainPODPNN');
    else
        save(['data/SVDResult_Ns', num2str(size(G, 1)),'_m',num2str(nModes), '_tol',num2str(tol),'_',num2str(freqSample(1)),'to',num2str(freqSample(length(freqSample))),'Hz'],'H','S','G', 'trainPODPNN');
    end
    clear DynamicSample
end

clear C Cpre Ccondpre 

%====================================================================================
% Online Stage
%====================================================================================
H=H(1:lenEM,:);
[nCond,~]=size(CondFactorOut);
DynamicTotal=cell(nCond);

stiffRed=H'*K(1:lenEM,1:lenEM)*H;
ResidRed=H'*Resid(1:lenEM);

stiffRedP=parallel.pool.Constant(stiffRed);
ResidRedP=parallel.pool.Constant(ResidRed);
HP=parallel.pool.Constant(H);
clear stiffRed
clear ResidRed


NmechBodies=ProblemData.NmechBodies;
for i=1:nCond
    Dynamic = zeros(nTotal,length(freqOut));
    C=CondFactorOut(i,1)*Ccond{1};
    for nBody=2:NmechBodies
    C=C+CondFactorOut(i,nBody)*Ccond{nBody};
    end
    dampRed=H'*C(1:lenEM,1:lenEM)*H-1i*H'*CReg(1:lenEM,1:lenEM)*H;
    dampRedP=parallel.pool.Constant(dampRed);
    clear dampRed
    parfor j=1:length(freqOut)
        
        % Define angular frequency (rad/s) from frequency (Hz)
        freq=freqOut(j);
        omega       = freq*2*pi;
        % Compute the K, C and M weights
        w0         = 1;
        w1         = complex(0,1)*omega;
        w2         = -omega^2;

        
        
        
        ReducedMatrix=w1*dampRedP.Value+stiffRedP.Value;
        Ared=ReducedMatrix\-ResidRedP.Value;
        AredWhole=HP.Value*Ared;
        %          AFull=linearSystemSolverEM(Unknown,ProblemData,w0,w1,CP.Value,CRegP.Value,KP.Value,ResidP.Value,CpreP.Value,CpreRegP.Value);
        % %         %DynamicEMFull(:,j)=AFull;
        %          ErrorVector(j)=norm(AredWhole-AFull(1:lenEM))/norm(AFull(1:lenEM));
        
        if Options.SplitMech==1
        [U] = linearSystemSolverMech(Unknown,AredWhole,w0,w1,w2,...
            MP.Value,KP.Value,ResidP.Value,dampRatio,NmechBodies);
        else
           [U] = linearSystemSolverMechNS(Unknown,AredWhole,w0,w1,w2,...
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
            disp(['------------ Conducting factor choice ',num2str(i),': ',num2str(j),'/',num2str(length(freqOut)), ' frequencies solved. ------------'])
            
    end
    DynamicTotal{i}=Dynamic;
end
% Rearrange solution into matrix
Dynamic=zeros(nTotal,nCond*length(freqOut));
for i=1:nCond
    Dynamic(:,(i-1)*length(freqOut)+1:i*length(freqOut))=DynamicTotal{i};
end

clear H U Ared AredWhole dampRed stiffRed ResidRed
clear M K Resid C COVC C77K C4K Cpre CpreReg CReg Cpre77K Cpre4K CpreOVC


