function [mass,damp,damppre,dampReg,damppreReg,stiff,Resid,Kd,Cd,CdReg,Md,nmst]=elementLoop(Static,StaticCurrent,Mesh,Basis,Quadrature,Unknown,UnknownCurrent,UnknownStatic,ProblemData,X,coupling,freqSweep)
%=========================================================================
% Extract the data from the data structure
%=========================================================================
% Extract the relevant data from the Mesh structure
intma           = Mesh.intma;
mat             = Mesh.mat;
nelem           = Mesh.Nelements;
coord           = Mesh.Coordinates;
subFlag         = Mesh.subFlag;
cond            = Mesh.face.cond;

% Extract relevant data from Unknown structure
nunktMech       = Unknown.Mech.nunkt;
ndirEM          = Unknown.EM.npec;
ndirM           = Unknown.Mech.npec;
nUnk            = Unknown.system.nunkt;
nSparse         = Unknown.system.nSparse;

% Extract relevant data from ProblemData structure
order           = ProblemData.jb.order;
sigma           = ProblemData.matr.sigma;
omega           = ProblemData.matr.omega;
muz             = ProblemData.matr.muz;
mu              = ProblemData.matr.mu;
regTerm         = ProblemData.matr.regTerm;
esizet          = ProblemData.esizet;
esizeH1         = ProblemData.esizeH1;
gorder          = ProblemData.jb.gorder;
probstatic      = ProblemData.probstatic;


%=========================================================================
% Begin initialisations
%=========================================================================

% Initialise arrays ready for assembly
Resid           = zeros(nUnk,1);
I               = zeros(nSparse,1);
J               = zeros(nSparse,1);
Kval            = zeros(nSparse,1);
Cval            = zeros(nSparse,1);
CvalReg          = zeros(nSparse,1);
Cvalpre         = zeros(nSparse,1);
CvalpreReg         = zeros(nSparse,1);
Mval            = [];
Kdval              =[];
Cdval              = [];
CdvalReg              = [];
Mdval              = [];
Idir             = [];
Jdir              = [];
% Initialise variables
count           = 0;
nmst            = 0;
nmstDir         = 0;

StressIntegral=zeros(3,1);
FlagMechVolume=Unknown.FlagMechVolume;
%=========================================================================
% Loop through elements
%=========================================================================
for i=1:nelem
    C1=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    CReg=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    Cpre=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    CpreReg=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    
    FlagVolume=FlagMechVolume(i);
    % Display the number of elements completed at intervals in the run
    count=count+1;
    if count==1000
        count=0;
        disp(['number of elements complete, nelem = ',num2str(i)]);
    end
    
    %---------------------------------------------------------------------
    % Get coordinates for this element based on vertices and blending
    %---------------------------------------------------------------------    
    % Coordinates of vertex
    xy = coord(intma(i,1:4),1:3);  
    
    % Blending function corrections (or not)
    if gorder>1
        mycoord=Mesh.mycoord(:,:,i);
        lec=Mesh.lec(:,:,:,i);
        lfc=Mesh.lfc(:,:,:,i);
        flag=Mesh.flagBlend(i);
    elseif gorder==1
        mycoord=Mesh.mycoord(:,:,i);  
        lec = zeros(6,gorder,3);
        lfc = zeros(4,gorder*(gorder-1)/2,3);
        flag=1;
    else
        gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;   
        mycoord = zeros(gesizet,3);
        mycoord(1:4,1:3)=xy;
        lec = zeros(6,gorder,3);
        lfc = zeros(4,gorder*(gorder-1)/2,3);
        flag=0;
    end

    %---------------------------------------------------------------------
    % Extract element type and material number from mesh structure
    %---------------------------------------------------------------------
    eltype=Mesh.eltype(i);
    material=mat(i);
    % Extract relative permeability for this element
    mue=mu(mat(i));
    
    %---------------------------------------------------------------------
    % Determine the unknown numbering of the elemental matrices
    %---------------------------------------------------------------------
    lunkv=Unknown.system.unknowns(:,i);
    lunkvEM=Unknown.EM.unknowns(:,i);
    lunkvM=Unknown.Mech.unknownsD(:,i);
    % Use static numbering to compute ADC
    lunkvEMStatic=UnknownStatic.EM.unknowns(:,i);
    lunkvEMCurrent=UnknownCurrent.EM.unknowns(:,i);

    
    %---------------------------------------------------------------------
    % Determine the elemental solution vector
    %---------------------------------------------------------------------
    Adc=zeros(esizet,1);
    A0=zeros(esizet,1);
    A=zeros(esizet,1);
    u=zeros(3*esizeH1,1);
    nunktEM=Unknown.EM.nunkt;
    
    % Extract the static solution (keep Dirichlet in case storage issue)
    Adc(lunkvEMStatic>0)=Static.sol(lunkvEMStatic(lunkvEMStatic>0),1);
    A0(lunkvEMCurrent>0)=StaticCurrent.sol(lunkvEMCurrent(lunkvEMCurrent>0),1);
    A(lunkvEM>0)=X(lunkvEM(lunkvEM>0),1);


    % Extract Mech solution (keep dirichlet in case storage issue)
    u(lunkvM>0)=X(lunkvM(lunkvM>0),1);  

    
    % Store the dynamic solution in the array
    Xdir  = [A;u];
    
    %---------------------------------------------------------------------
    % Determine the Elemental matrices
    %---------------------------------------------------------------------
    [M,C,K,Res,StressIntegral]=GaussIntegStatic(Basis,Quadrature,ProblemData,...
        flag,mue,Adc,A0,A,u,mycoord,eltype,probstatic,material,subFlag(i),...
        i,lec,lfc,xy,cond,coupling,StressIntegral,FlagVolume);

    
    %-------------------------------------------------------------------------
    % Define kappa coefficient defined by regularization term and conductivity
    %-------------------------------------------------------------------------
    
    if subFlag(i)==2                                   % Conducting region
        if probstatic==0                            % Dynamic problem
            kappalow=sigma(mat(i))*muz;
            kappalowReg=0;
        else                     % Static problem
            kappalow=0; 
            kappalowReg=regTerm;
        end
    else   
          kappalow=0;
          kappalowReg=regTerm;
    end
    
    % Redefine regularisation for solver option 1
    if ProblemData.sol.regopt ==1
        if subFlag(i)~=2
            kappa=0;
            kappaReg=0;
        else
            kappa=kappalow;
            kappaReg=kappalowReg;
        end
    else
        kappa=kappalow;
        kappaReg=kappalowReg;
    end

    
    % Regularization
    C1(1:6,1:6)=kappalow*C(1:6,1:6);
    CReg(1:6,1:6)=kappalowReg*C(1:6,1:6);
    Cpre(1:6,1:6)=C1(1:6,1:6);
    CpreReg(1:6,1:6)=CReg(1:6,1:6);

     

    % Extract number of interior EM basis
    nintbas=Unknown.EM.nintbas;
    
    
    %  Add regularisation to higher order block
    for p=1:esizet
        if p<= 6
            C1(p,7:esizet)=kappa*C(p,7:esizet);
            CReg(p,7:esizet)=kappaReg*C(p,7:esizet);
        else
            C1(p,1:esizet)=kappa*C(p,1:esizet);
            CReg(p,1:esizet)=kappaReg*C(p,1:esizet);
            Cpre(p,1:esizet)=abs(kappa)*C(p,1:esizet);
            CpreReg(p,1:esizet)=abs(kappaReg)*C(p,1:esizet);
        end
    end
    
    if ProblemData.sol.regopt==1
        % high order edges (gradients)
        Cpre(7:6*(order+1),7:6*(order+1))=abs(kappa)*C(7:6*(order+1),7:6*(order+1));
        CpreReg(7:6*(order+1),7:6*(order+1))=abs(kappaReg)*C(7:6*(order+1),7:6*(order+1));
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        Cpre(st:en,st:en)=abs(kappa)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappaReg)*C(st:en,st:en);
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        Cpre(st:en,st:en)=abs(kappalow)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        Cpre(st:en,st:en)=abs(kappa)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappaReg)*C(st:en,st:en);
        
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        Cpre(st:en,st:en)=abs(kappalow)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
    else
        % high order edges (gradients)
        Cpre(7:6*(order+1),7:6*(order+1))=abs(kappalow)*C(7:6*(order+1),7:6*(order+1));
        CpreReg(7:6*(order+1),7:6*(order+1))=abs(kappalowReg)*C(7:6*(order+1),7:6*(order+1));
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        Cpre(st:en,st:en)=abs(kappalow)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        Cpre(st:en,st:en)=abs(kappalow)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        Cpre(st:en,st:en)=abs(kappalow)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        Cpre(st:en,st:en)=abs(kappalow)*C(st:en,st:en);
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
    end
   
    
    % Complete damping matrix (no regularisation for mechanic block)
    C1(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1)=C(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1);
    Cpre(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1)=C(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1);
    C1(1:esizet,esizet+1:esizet+3*esizeH1)=C(1:esizet,esizet+1:esizet+3*esizeH1);
    Cpre(1:esizet,esizet+1:esizet+3*esizeH1)=C(1:esizet,esizet+1:esizet+3*esizeH1);
    %---------------------------------------------------------------------
    % Linear system assembly in vector format
    %---------------------------------------------------------------------
    [Mval,Cval,Cvalpre,CvalReg,CvalpreReg,Kval,Resid,I,J,Idir,Jdir,nmst,nmstDir,Kdval,Cdval,CdvalReg,Mdval]=matrixAssembly(K,C1,Cpre,CReg,CpreReg,M,Res,Xdir,Mval,Cval,Cvalpre,CvalReg,CvalpreReg,Kval,Resid,I,J,Idir,Jdir,lunkv,probstatic,esizet,esizeH1,nmst,nmstDir,omega,Kdval,Cdval,CdvalReg,Mdval,freqSweep);
    
end


%=========================================================================
% Construct the Linear system in sparse format
%=========================================================================
% Build stiffness and damping matrices from I,J,Kval and Cval
stiff  = sparse(I(1:nmst),J(1:nmst),Kval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damp  = sparse(I(1:nmst),J(1:nmst),Cval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
dampReg  = sparse(I(1:nmst),J(1:nmst),CvalReg(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
mass  = sparse(I(1:nmst),J(1:nmst),Mval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damppre=sparse(I(1:nmst),J(1:nmst),Cvalpre(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damppreReg=sparse(I(1:nmst),J(1:nmst),CvalpreReg(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
Kd=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),Kdval(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
Cd=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),Cdval(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
CdReg=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),CdvalReg(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
Md=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),Mdval(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
end


