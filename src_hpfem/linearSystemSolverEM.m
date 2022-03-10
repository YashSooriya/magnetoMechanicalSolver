function [U]=linearSystemSolverEM(Unknown,ProblemData,w0,w1,w2,...
    M,C,CReg,K,Resid,Cpre,CpreReg,Options)

%=========================================================================
% Initialisations
%=========================================================================
lenEM=Unknown.EM.nunkt;
ndirEM=Unknown.EM.npec;


%=========================================================================
% System assembly
%=========================================================================
wReg=(-1i)*w1;

% Construct the Tangent stiffness matrix
TSM = w2*M+w1*C +wReg*CReg+ w0*K;

clear C
clear CReg
% Construct matrix for preconditioner
TSMpre = w2*M+w1*Cpre+wReg*CpreReg+w0*K;

clear Cpre
clear K
clear CpreReg

%=========================================================================
% Re-structure the matrix blocks to solve as fixed point
%=========================================================================
% extract all the blocks: AA, AU, UA, UU from TSM

TSM_AA=TSM(1:lenEM,1:lenEM);
TSM_AApre=TSMpre(1:lenEM,1:lenEM);
clear TSMpre

clear TSM
clear M
tic

%=========================================================================
% Preconditioning for iterative solver
%=========================================================================
if ProblemData.sol.regopt~=3      % Iterative solvers
    
    % Set up data for the preconditioner
    [ Xzz,XEgEg,XFgFg, XFF, XIgIg, XII, nz, nheg, nhfg, nhf,nhig, nhi] = extract_em(Unknown,TSM_AApre,TSM_AA,ProblemData);
    tol=ProblemData.TOL_GMRES;
    maxit=1500;
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
clear TSM_AApre
 
    Resid_EM=Resid(1:lenEM);
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
    
    elseif ProblemData.sol.regopt==2
    % Solve linear system with direct LU solve for application of the
    % action of the inverse of the gradient blocks
    
    I=[1:lenEM];
    J=[1:lenEM];
    X1=ones(1,lenEM);
    preconlm=sparse(I,J,X1,lenEM,lenEM);
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
end


    
%-------------------------------------------------------------------------
% Solve the linear EM system with iterative or direct solver
%-------------------------------------------------------------------------

if ProblemData.sol.regopt ==1
    
    % Solve linear system with preconditioned GMRES    
    [U_AA,flag,~,~,resvec] = gmres(TSM_AA,-Resid_EM,10,tol,maxit,preconlm,@(x)preconditioner(x,arg));
    
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
        
    % Solve linear system with preconditioned GMRES
    [U_AA,flag,relres,iter,resvec] = gmres(TSM_AA,-Resid_EM,10,tol,maxit,preconlm,@(x)preconditioner_dir(x,arg));
    
    
    % Display flag if there is a problem with GMRES
    if flag~=0
        disp('problem with gmres')
        flag
        
    else
        disp('System solved succesfully');
        disp(['Number of iterations',num2str(length(resvec))]);
        
    end
    
else
   
    % Direct Solver
     U_AA= TSM_AA\-Resid_EM;
end


% order=ProblemData.jb.order+1;
% fileName=sprintf('IterationsGMRESToyToleranceMinus5_1000Hz_p%d',order);
% save(fileName,'resvec');
toc
disp('for solving the linear EM system with iterative or direct solver.')
U=[U_AA;zeros(ndirEM,1)];


end