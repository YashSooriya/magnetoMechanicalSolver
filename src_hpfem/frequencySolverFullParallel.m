function [Dynamic, Toc, toccount]=frequencySolverFullParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,Damp,dampRatio,CondFactorOut,CondFactorChoice,Ncores, ticInit, Toc, toccount)
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

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

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

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%-------------------------------------------------------------------------
% Frequency sweep
%-------------------------------------------------------------------------

%---------------------------------------------------------------------
% Linear system solver
%---------------------------------------------------------------------
% Include Dirichlet values in solution
if Options.Non0Dir==1 && freqSweep==0
    % Compute BC coefficents if they are non-zero and not changing with
    % frequency
    % Initial guess of the solution field
    ProblemData.probFlag=1;
    [AAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.EM,ProblemData,2*pi*freqOut);
    % Define the problem flag
    ProblemData.probFlag=2;
    
    % Initial guess for the mechanical problem
    [UAC] = initialGuess(Mesh,Basis,Quadrature,Unknown.Mech,ProblemData,2*pi*freqOut);
    X=[AAC;UAC];
else
    X=zeros(nTotal,1);
end

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

%=========================================================================
% System assembly
%=========================================================================
fileName=sprintf('MatricesDynamic_%d',ProblemData.order);
fileName=strcat(fileName,ProblemData.jb.job);
if Assemble==1
    disp('Constructing the Stiffness and Forcing Matrices')
    % Assemble the system components into the tangent matrix and Residual
    [M,Ccond,Ccondpre,CReg,CpreReg,K,Resid_fixed,Kd,Cdcond,CdReg,Md,~]=elementLoop(Static,StaticCurrent,Mesh,Basis,Quadrature,Unknown,UnknownCurrent,UnknownStatic,ProblemData,X,couple,freqSweep,SourceMapping);
    %save(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','Kd','Cd','CdReg','Md','-v7.3');
else
    load(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','Kd','Cd','CdReg','Md', 'Ccond', 'Ccondpre', 'Resid');
    %Resid=Res;
end

% Legacy code note that X has only 1 column and so this will not work.
% if freqSweep==1
%     dirDOF=lenEM+[(1:ndirEM) ndirEM+lenM+(1:ndirMech)];
%     % Add Dirichlet contributions to residual vector (the vector X contains the weights);
%     Resid=Resid+Md*X(dirDOF,2)+Cd*X(dirDOF,2)+(Kd+wReg*CdReg)*X(dirDOF,1);
% end


% Get number of different conductivities to solve for
[nCond,~]=size(CondFactorOut);
% Initialize the solution matrix
DynamicTotal=cell(nCond);

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;    
parToc = zeros(1,length(freqOut));

% Start parallel pool
parpool(Ncores)
% Copy matrices to workers for parallel efficiency
MP=parallel.pool.Constant(M);
KP=parallel.pool.Constant(K);
ResidP=parallel.pool.Constant(Resid_fixed);
CRegP=parallel.pool.Constant(CReg);
CpreRegP=parallel.pool.Constant(CpreReg);

Toc(toccount) = toc(ticInit);
toccount = toccount + 1;

for i=1:nCond
    Dynamic=zeros(nTotal,length(freqOut));
    C=CondFactorOut(i,1)*Ccond{1};
    Cpre=CondFactorOut(i,1)*Ccondpre{1};
    for nBody=2:NmechBodies
    C=C+CondFactorOut(i,nBody)*Ccond{nBody};
    Cpre=Cpre+CondFactorOut(i,nBody)*Ccondpre{nBody};
    end
    %We want to Cd which are the contributions to the damping matrix
    %relating to Dircihlet DOFs. These need to be built in the same way as
    %C so we want to allow for each mechanical body to have a different
    %conductivity factor.
    Cd=CondFactorOut(i,1)*Cdcond{1};
    for nBody=2:NmechBodies
    Cd=Cd+CondFactorOut(i,nBody)*Cdcond{nBody};
    end

    CP=parallel.pool.Constant(C);
    CpreP=parallel.pool.Constant(Cpre);

    Toc(toccount) = toc(ticInit);
    toccount = toccount + 1;

    parfor j=1:length(freqOut)
        % Define angular frequency (rad/s) from frequency (Hz)
        freq=freqOut(j);
        omega       = freq*2*pi;
        % Compute the K, C and M weights
        w0         = 1;
        w1         = complex(0,1)*omega;
        w2         = -omega^2;
        wReg       = (-1i)*w1;
        

        
        
        
        
        if Options.Non0Dir==1 && freqSweep==1
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

        dirDOF=lenEM+[(1:ndirEM) ndirEM+lenM+(1:ndirMech)];
        Resid=Resid_fixed+w2*Md*X(dirDOF)+w1*Cd*X(dirDOF)+(w0*Kd+wReg*CdReg)*X(dirDOF);

        if Options.Non0Dir==1 && freqSweep==1
            [U] = linearSystemSolver(Unknown,ProblemData,w0,w1,w2,...
                MP.Value,CP.Value,CRegP.Value,KP.Value,Resid,CpreP.Value,CpreRegP.Value,Options,dampRatio);
        else
            [U] = linearSystemSolverNS(Unknown,ProblemData,w0,w1,w2,...
                MP.Value,CP.Value,CRegP.Value,KP.Value,Resid,CpreP.Value,CpreRegP.Value,dampRatio);
        end

        % Update solution vector (if Dirichlet values contained in X, otherwise X=0 and X=X+U=U);
        U=U+X;

        %---------------------------------------------------------------------
        % Store the solution at each frequency
        %---------------------------------------------------------------------
        
        Dynamic(:,j)= U;

    end

    DynamicTotal{i}=Dynamic;
end

freqToc = toc(ticInit) - Toc(toccount - 1);
freqTocMean = freqToc / length(freqOut);

for ii = 1:length(freqOut)
    Toc(toccount) = Toc(toccount - 1) + freqTocMean;
    toccount = toccount + 1;
end

% Rearrange solution into matrix
Dynamic=zeros(nTotal,nCond*length(freqOut));
for i=1:nCond
    Dynamic(:,(i-1)*length(freqOut)+1:i*length(freqOut))=DynamicTotal{i};
end




