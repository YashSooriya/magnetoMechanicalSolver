function [U,C,K,Res,CondNum]=linearSystemSolver2(Static,mesh,Basis,Quadrature,unknown,ProblemData,X,w0,w1,w2,coupling,freqSweep,Assemble)

%=========================================================================
% Initialisations
%=========================================================================
lenEM=unknown.EM.nunkt;
lenM=unknown.Mech.nunkt;
ndirEM=unknown.EM.npec;
ndirM=unknown.Mech.npec;

%=========================================================================
% System assembly
%=========================================================================
disp('Constructing the Stiffness and Forcing Matrices')
tic
% Assemble the system components into the tangent matrix and Residual
fileName=sprintf('MatricesMech_%d',ProblemData.order);
fileName=strcat(fileName,ProblemData.jb.job);
if Assemble==1
[M,C,~,CReg,~,K,Res,~,~,~,~,~]=elementLoop(Static,mesh,Basis,Quadrature,unknown,unknown,ProblemData,X,coupling,freqSweep);
save(fileName,'M','C','CReg','K','Res');
else
 load(fileName,'M','C','CReg','K','Res');  
end


% Construct the Tangent stiffness matrix
TSM = w2*M+w1*C +CReg+ w0*K;


% Solve just the mechanical part (EM solved previously)
nunkt=unknown.system.nunkt;
lenM=nunkt-lenEM;
TSM2=TSM(lenEM+1:end,lenEM+1:end);
Res2=Res(lenEM+1:end);


% % Compute the estimated condition number
%CondNum=condest(TSM);
%disp(['Condition Number of TSM, cond(K)=',num2str(CondNum)])
CondNum=0;
%=========================================================================
% Compute the increment vector U by solving the coupled linear system
%=========================================================================
% Solve the linear system

if ProblemData.sol.regopt~=3    % Iterative solvers
    
    % Set up data for the preconditioner
    [stiffVV,stiffEE,stiffFF,stiffII,nVV,nEE,nFF,nII] = extract_mech2(unknown,TSM2,ProblemData);
    tol=1e-7;
    maxit=50000;
    
    arg.XVV=stiffVV;
    arg.XEE=stiffEE;
    arg.XFF=stiffFF;
    arg.XII=stiffII;
    arg.nVV=nVV;
    arg.nEE=nEE;
    arg.nFF=nFF;
    arg.nII=nII;
    arg.nunkt=lenEM;
    
    
% %=========================================================================
% % Build preconditioner matrix for eigenspectrum analysis
% %=========================================================================
% PrecondSize=nVV+nEE+nFF+nII;
% PrecondMatrix=zeros(PrecondSize);
% PrecondMatrix(1:nVV,1:nVV)=stiffVV;
% PrecondMatrix(nVV+1:nVV+nEE,nVV+1:nVV+nEE)=stiffEE;
% PrecondMatrix(nVV+nEE+1:nVV+nEE+nFF,nVV+nEE+1:nVV+nEE+nFF)=stiffFF;
% PrecondMatrix(nVV+nEE+nFF+1:nVV+nEE+nFF+nII,nVV+nEE+nFF+1:nVV+nEE+nFF+nII)=stiffII;
% 
% eigenValues=eig(PrecondMatrix\TSM2);
% eigReal=real(eigenValues);
% eigImag=imag(eigenValues);
% TSM3=full(TSM2);
% eigenValues2=eig(TSM3);
% eigReal2=real(eigenValues2);
% eigImag2=imag(eigenValues2);
% figure()
% plot(eigReal,eigImag,'b*')
% figure
% plot(eigReal2,eigImag2,'b*')
%    aaaaa=1; 
    
end
if ProblemData.sol.regopt ==1
    
    % Solve linear system with iterative solve for application of the
    % action of the inverse of the gradient blocks
    I=[1:lenM];
    J=[1:lenM];
    X=ones(1,lenM);
    preconlm=sparse(I,J,X,lenM,lenM);
    
    % LU decomposition for mechanical blocks
    [Lm,Um,Pm,Qm] = lu(stiffVV);
    disp('completed LU decomposition')
    arg.Lm=Lm;
    arg.Um=Um;
    arg.Pm=Pm;
    arg.Qm=Qm;
    if nEE> 0
        [Lem,Uem,Pem,Qem] = lu(stiffEE);
        disp('completed LU decomposition')
        arg.Lem=Lem;
        arg.Uem=Uem;
        arg.Pem=Pem;
        arg.Qem=Qem;
    end
    
    if nFF > 0
        [Lfm,Ufm,Pfm,Qfm] = lu(stiffFF);
        disp('completed LU decomposition')
        arg.Lfm=Lfm;
        arg.Ufm=Ufm;
        arg.Pfm=Pfm;
        arg.Qfm=Qfm;
    end
    
    if nII > 0
        [Lim,Uim,Pim,Qim] = lu(stiffII);
        disp('completed LU decomposition')
        arg.Lim=Lim;
        arg.Uim=Uim;
        arg.Pim=Pim;
        arg.Qim=Qim;
    end
    
    % Solve the system with preconditioned GMRES
    [U,flag,relres,iter,resvec] = gmres(TSM2,-Res2,10,tol,maxit,preconlm,@(x)preconditioner2(x,arg));
    
    % Display flag if there is a problem with GMRES
    if flag~=0
        disp('problem with gmres')
        flag
        
    else
        disp('System solved succesfully');
        disp(['Number of iterations',num2str(length(resvec))]);
        
    end
    
elseif ProblemData.sol.regopt==2
    % Solve linear system with direct LU solve for application of the
    % action of the inverse of the gradient blocks
    
    I=[1:lenM];
    J=[1:lenM];
    X=ones(1,lenM);
    preconlm=sparse(I,J,X,lenM,lenM);
    
    % LU decomposition for mechanical blocks
    [Lm,Um,Pm,Qm] = lu(stiffVV);
    disp('completed LU decomposition')
    arg.Lm=Lm;
    arg.Um=Um;
    arg.Pm=Pm;
    arg.Qm=Qm;
    if nEE> 0
        [Lem,Uem,Pem,Qem] = lu(stiffEE);
        disp('completed LU decomposition')
        arg.Lem=Lem;
        arg.Uem=Uem;
        arg.Pem=Pem;
        arg.Qem=Qem;
    end
    
    if nFF > 0
        [Lfm,Ufm,Pfm,Qfm] = lu(stiffFF);
        disp('completed LU decomposition')
        arg.Lfm=Lfm;
        arg.Ufm=Ufm;
        arg.Pfm=Pfm;
        arg.Qfm=Qfm;
    end
    
    if nII > 0
        [Lim,Uim,Pim,Qim] = lu(stiffII);
        disp('completed LU decomposition')
        arg.Lim=Lim;
        arg.Uim=Uim;
        arg.Pim=Pim;
        arg.Qim=Qim;
    end
    
    % Solve the system with preconditioned GMRES
    [U,flag,relres,iter,resvec] = gmres(TSM2,-Res2,10,tol,maxit,preconlm,@(x)preconditioner_dir2(x,arg));
    
    if flag~=0
        disp('problem with gmres')
        flag
        
    else
        disp('System solved succesfully');
        disp(['Number of iterations',num2str(length(resvec))]);
        
    end
    
else
    
    % Direct Solve
    U=TSM2\-Res2;
    
end

%=========================================================================
% Save number of iterations required to converge
%=========================================================================
order=ProblemData.orderH1;
% fileName=sprintf('IterationsMechanicsStatic1838elem_p%d',order);
% save(fileName,'resvec');
%=========================================================================


% Complete the update with zeros for the electromagnetic part
U=[zeros(lenEM,1);U];

toc

% Extract the unknown numbering data
ndirEM   = unknown.EM.npec;
lenMech  = unknown.Mech.nunkt;
ndirM = unknown.Mech.npec;


% Include the Dirichlet Dof into the update vector U, where Dir=0;
U=[U(1:lenEM);zeros(ndirEM,1);U(lenEM+1:lenMech+lenEM);zeros(ndirM,1)];

end