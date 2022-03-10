function [U,CondNum,unknown]=linearSystemSolver1NS(Static,StaticCurrent,mesh,Basis,Quadrature,unknown,ProblemData,X,w0,w1,w2,coupling,freqSweep,Assemble,staticMechanics,SourceMapping)

%=========================================================================
% Initialisations
%=========================================================================
lenEM=unknown.EM.nunkt;
coupling=Options.Couple;
staticMechanics=Options.StaticMechanics;
freqSweep=Options.freqSweep;
Assemble=Options.Assemble;
SourceMapping=Options.SourceMapping;


%=========================================================================
% System assembly
%=========================================================================
disp('Constructing the Stiffness and Forcing Matrices')

% Assemble the system components into the tangent matrix and Residual
tic
fileName=sprintf('MatricesEM_%d',ProblemData.order);
fileName=strcat(fileName,ProblemData.jb.job);
if Assemble==1
  [M,COVC,C77K,C4K,CpreOVC,Cpre77K,Cpre4K,CReg,CpreReg,K,Res,~,~,~,~,nmst]=elementLoop(Static,StaticCurrent,mesh,Basis,Quadrature,unknown,unknown,unknown,ProblemData,X,coupling,freqSweep,SourceMapping);
    %save(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','nmst','-v7.3');
else
    load(fileName,'M','C','Cpre','CReg','CpreReg','K','Res','nmst');
end
% Overwrite the estimated size of sparse matrix vectors with exact value
% Not overallocating or underallocating for next solver iterations
unknown.system.nSparse=nmst;



toc
TimeAssembly=toc

% Construct the Tangent stiffness matrix
TSM = w2*M+w1*(COVC+C77K+C4K) +CReg+ w0*K;
TSMpre = w2*M+w1*CpreOVC+w1*Cpre77K+w1*Cpre4K+CpreReg+w0*K;



% Extract EM matrix
TSM2=TSM(1:lenEM,1:lenEM);
Res2=Res(1:lenEM);
TSMpre2=TSMpre(1:lenEM,1:lenEM);
if staticMechanics==1
    % Extract mechanical matrix
    nunkt=unknown.system.nunkt;
    lenM=nunkt-lenEM;
    TSM_UU=TSM(lenEM+1:end,lenEM+1:end);
    Res_U=Res(lenEM+1:end);
end
clear TSM
clear TSMpre

CondNum=0;
%=========================================================================
% Compute the increment vector U by solving the coupled linear system
%=========================================================================
% Solve the linear system
tic

if ProblemData.sol.regopt~=3      % Iterative solvers
    
    % Set up data for the preconditioner
    [ Xzz,XEgEg,XFgFg, XFF, XIgIg, XII, nz, nheg, nhfg, nhf,nhig, nhi] = extract_mech(unknown,TSMpre2,TSM2,ProblemData);
    tol=1e-7;
    maxit=3000;
    arg.Xzz=Xzz;
    arg.nz=nz;
    arg.XEgEg=XEgEg;
    arg.nheg=nheg;
    arg.XFgFg=XFgFg;
    arg.nhfg=nhfg;
    arg.XFF=XFF;
    arg.nhf=nhf;
    arg.XIgIg=XIgIg;
    arg.nhig=nhig;
    arg.XII=XII;
    arg.nhi=nhi;
    arg.nunkt=lenEM;
end

if ProblemData.sol.regopt ==1
    % Solve linear system with iterative solve for application of the
    % action of the inverse of the gradient blocks
    
    I=[1:lenEM];
    J=[1:lenEM];
    X=ones(1,lenEM);
    preconlm=sparse(I,J,X,lenEM,lenEM);
    
    [L,U,P,Q] = lu(Xzz);
    disp('completed LU decomposition')
    arg.L=L;
    arg.U=U;
    arg.P=P;
    arg.Q=Q;
    if nhf > 0
        [Lf,Uf,Pf,Qf] = lu(XFF);
        disp('completed LU decomposition')
        arg.Lf=Lf;
        arg.Uf=Uf;
        arg.Pf=Pf;
        arg.Qf=Qf;
    end
    
    if nhi > 0
        [Li,Ui,Pi,Qi] = lu(XII);
        disp('completed LU decomposition')
        arg.Li=Li;
        arg.Ui=Ui;
        arg.Pi=Pi;
        
        arg.Qi=Qi;
    end
    if nhig>0
        [Lig,Uig,Pig,Qig] = lu(XIgIg);
        disp('completed LU decomposition')
        arg.Lig=Lig;
        arg.Uig=Uig;
        arg.Pig=Pig;
        arg.Qig=Qig;
    end
    
    % Solve linear system with preconditioned GMRES
    [U,flag,~,~,resvec] = gmres(TSM2,-Res2,10,tol,maxit,preconlm,@(x)preconditioner2(x,arg));
    
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
    
    I=[1:lenEM];
    J=[1:lenEM];
    X=ones(1,lenEM);
    preconlm=sparse(I,J,X,lenEM,lenEM);
    [L,U,P,Q] = lu(Xzz);
    disp('completed LU decomposition')
    arg.L=L;
    arg.U=U;
    arg.P=P;
    arg.Q=Q;
    if nheg > 0
        [Leg,Ueg,Peg,Qeg] = lu(XEgEg);
        disp('completed LU decomposition')
        arg.Leg=Leg;
        arg.Ueg=Ueg;
        arg.Peg=Peg;
        arg.Qeg=Qeg;
    end
    if nhfg > 0
        [Lfg,Ufg,Pfg,Qfg] = lu(XFgFg);
        disp('completed LU decomposition')
        arg.Lfg=Lfg;
        arg.Ufg=Ufg;
        arg.Pfg=Pfg;
        arg.Qfg=Qfg;
    end
    if nhf > 0
        [Lf,Uf,Pf,Qf] = lu(XFF);
        disp('completed LU decomposition')
        arg.Lf=Lf;
        arg.Uf=Uf;
        arg.Pf=Pf;
        arg.Qf=Qf;
    end
    if nhig>0
        [Lig,Uig,Pig,Qig] = lu(XIgIg);
        disp('completed LU decomposition')
        arg.Lig=Lig;
        arg.Uig=Uig;
        arg.Pig=Pig;
        arg.Qig=Qig;
    end
    
    if nhi > 0
        [Li,Ui,Pi,Qi] = lu(XII);
        disp('completed LU decomposition')
        arg.Li=Li;
        arg.Ui=Ui;
        arg.Pi=Pi;
        arg.Qi=Qi;
    end
    
    % Solve linear system with preconditioned GMRES
    [U,flag,relres,iter,resvec] = gmres(TSM2,-Res2,10,tol,maxit,preconlm,@(x)preconditioner_dir(x,arg));
    
    
    % Display flag if there is a problem with GMRES
    if flag~=0
        disp('problem with gmres')
        flag
        
    else
        disp('System solved succesfully');
        disp(['Number of iterations',num2str(length(resvec))]);
        
    end
    
else
    
    % Direct Solve
    U= TSM2\-Res2;
    
end

% Set mechanical field solution to zero (will update afterwards)
nunkt=unknown.system.nunkt;
U(lenEM+1:nunkt)=zeros(nunkt-lenEM,1);
toc


% Extract the unknown numbering data
ndirEM   = unknown.EM.npec;
lenMech  = unknown.Mech.nunkt;
ndirMech = unknown.Mech.npec;


% Include the Dirichlet Dof into the update vector U, where Dir=0;
U=[U(1:lenEM);zeros(ndirEM,1);U(lenEM+1:lenMech+lenEM);zeros(ndirMech,1)];


%============================================================================================================
% Solve the mechanical DC problem
%============================================================================================================
if staticMechanics==1
    if ProblemData.sol.regopt~=3    % Iterative solvers
        
        % Set up data for the preconditioner
        [stiffVV,stiffEE,stiffFF,stiffII,nVV,nEE,nFF,nII] = extract_mech2(unknown,TSM_UU,ProblemData);
        tol=1e-6;
        maxit=5000;
        
        arg.XVV=stiffVV;
        arg.XEE=stiffEE;
        arg.XFF=stiffFF;
        arg.XII=stiffII;
        arg.nVV=nVV;
        arg.nEE=nEE;
        arg.nFF=nFF;
        arg.nII=nII;
        arg.nunkt=lenEM;
        
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
        [U2,flag,relres,iter,resvec] = gmres(TSM_UU,-Res_U,10,tol,maxit,preconlm,@(x)preconditioner2(x,arg));
        
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
        [U2,flag,relres,iter,resvec] = gmres(TSM_UU,-Res_U,10,tol,maxit,preconlm,@(x)preconditioner_dir2(x,arg));
        
        if flag~=0
            disp('problem with gmres')
            flag
            
        else
            disp('System solved succesfully');
            disp(['Number of iterations',num2str(length(resvec))]);
            
        end
        
    else
        
        % Direct Solve
        U2=TSM_UU\-Res_U;
        
    end
    
    
    
    % Complete the update with zeros for the electromagnetic part
    U2=[zeros(lenEM,1);U2];
    
    toc
    
    % Extract the unknown numbering data
    ndirEM   = unknown.EM.npec;
    lenMech  = unknown.Mech.nunkt;
    ndirM = unknown.Mech.npec;
    
    
    % Include the Dirichlet Dof into the update vector U, where Dir=0;
    U=[U(1:lenEM);zeros(ndirEM,1);U2(lenEM+1:lenMech+lenEM);zeros(ndirM,1)];
    
end
