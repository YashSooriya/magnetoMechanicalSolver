function [Dynamic]=frequencySolverFullSerial(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,Damp,dampRatio,CondFactorOut,CondFactorChoice)
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


%=========================================================================
% System assembly
%=========================================================================
fileName=sprintf('MatricesDynamic_%d',ProblemData.order);
fileName=strcat(fileName,ProblemData.jb.job);
if Assemble==1
    disp('Constructing the Stiffness and Forcing Matrices')
    % Assemble the system components into the tangent matrix and Residual
    [M,Ccond,Ccondpre,CReg,CpreReg,K,Resid,Kd,Cd,CdReg,Md,~]=elementLoop(Static,StaticCurrent,Mesh,Basis,Quadrature,Unknown,UnknownCurrent,UnknownStatic,ProblemData,X,couple,freqSweep,SourceMapping);
    %save(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','Kd','Cd','CdReg','Md','-v7.3');
else
    load(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','Kd','Cd','CdReg','Md', 'Ccond', 'Ccondpre', 'Resid');
    %Resid=Res;
end

% Legacy code note that X has only 1 column and so this will not work.
% if freqSweep==1
%     dirDOF=lenEM+[(1:ndirEM) ndirEM+lenM+(1:ndirMech)];
%     % Add Dirichlet contributions to residual vector (the vector X contains the weights);
% %     Resid=Resid+Md*X(dirDOF,2)+Cd*X(dirDOF,2)+(Kd+wReg*CdReg)*X(dirDOF,1);
%     Resid=Resid+Md*X(dirDOF)+Cd*X(dirDOF)+(Kd+wReg*CdReg)*X(dirDOF);
% end


% Get number of different conductivities to solve for
[nCond,~]=size(CondFactorOut);
% Initialize the solution matrix
Dynamic = zeros(nTotal,nCond*length(freqOut));
% Initialize counter
count=0;
for i=1:nCond
    C=CondFactorOut(i,1)*Ccond{1};
    Cpre=CondFactorOut(i,1)*Ccondpre{1};
    for nBody=2:NmechBodies
    C=C+CondFactorOut(i,nBody)*Ccond{nBody};
    Cpre=Cpre+CondFactorOut(i,nBody)*Ccondpre{nBody};
    end
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
            ProblemData.probFlag=1;
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
            [U] = linearSystemSolver(Unknown,ProblemData,w0,w1,w2,...
                M,C,CReg,K,Resid,Cpre,CpreReg,Options,dampRatio);
        else
            [U] = linearSystemSolverNS(Unknown,ProblemData,w0,w1,w2,...
                M,C,CReg,K,Resid,Cpre,CpreReg,dampRatio);
        end
        % Update solution vector (if Dirichlet values contained in X, otherwise X=0 and X=X+U=U);
        U=U+X;
        %---------------------------------------------------------------------
        % Store the solution at each frequency
        %---------------------------------------------------------------------
        
        Dynamic(:,count)= U;
        
    end
end



