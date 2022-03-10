function [U]=linearSystemSolver3(Static,StaticCurrent,mesh,Basis,Quadrature,unknown,ProblemData,X,w0,w1,w2,Options)

%=========================================================================
% Initialisations
%=========================================================================
lenEM=unknown.EM.nunkt;
coupling=Options.Couple;
freqSweep=Options.freqSweep;
Assemble=Options.Assemble;
SourceMapping=Options.SourceMapping;


%=========================================================================
% System assembly
%=========================================================================
disp('Constructing the Stiffness and Forcing Matrices')

% Assemble the system components into the tangent matrix and Residual
tic
fileName=sprintf('MatricesEMCurrent_%d',ProblemData.order);
fileName=strcat(fileName,ProblemData.jb.job);
if Assemble==1
  [~,~,~,CReg,CpreReg,K,Res,~,~,~,~,nmst]=elementLoop(Static,StaticCurrent,mesh,Basis,Quadrature,unknown,unknown,unknown,ProblemData,X,coupling,freqSweep,SourceMapping);
%save(fileName,'CReg','CpreReg','K','Res','nmst','-v7.3');
else
load(fileName,'CReg','CpreReg','K','Res','nmst');
end
% Overwrite the estimated size of sparse matrix vectors with exact value
% Not overallocating or underallocating for next solver iterations
unknown.system.nSparse=nmst;



toc
TimeAssembly=toc

% Construct the Tangent stiffness matrix
TSM = 10*CReg+ w0*K;  
TSMpre =10*CpreReg+w0*K;


% Solve just the EM part (Mechanics is solved afterwards)
TSM2=TSM(1:lenEM,1:lenEM);
Res2=Res(1:lenEM);
TSMpre2=TSMpre(1:lenEM,1:lenEM);
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
    tol=1e-9;
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
    [U,flag,relres,iter,resvec] = gmres(TSM2,-Res2,50,tol,maxit,preconlm,@(x)preconditioner_dir(x,arg));
    %save('residualCoil','resvec');
    
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

end