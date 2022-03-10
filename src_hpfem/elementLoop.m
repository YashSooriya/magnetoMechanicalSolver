function [mass,damp,damppre,dampReg,damppreReg,stiff,Resid,Kd,Cd,CdReg,Md,nmst]=elementLoop(Static,StaticCurrent,Mesh,Basis,Quadrature,Unknown,UnknownCurrent,UnknownStatic,ProblemData,X,coupling,freqSweep,SourceMapping)
%=========================================================================
% Extract the data from the data structure
%=========================================================================
% Extract the relevant data from the Mesh structure
intma           = Mesh.intma;
mat             = Mesh.mat;
nelem           = Mesh.Nelements;
coord           = Mesh.Coordinates;
matc            = Mesh.matc;
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
NmechBodies     = ProblemData.NmechBodies;
matMech         =ProblemData.matMech;

%=========================================================================
% Begin initialisations
%=========================================================================

% Initialise arrays ready for assembly
Resid           = zeros(nUnk,1);
I               = zeros(nSparse,1);
J               = zeros(nSparse,1);
Kval            = zeros(nSparse,1);
Cval            = cell(NmechBodies,1);
Cvalpre         = cell(NmechBodies,1);
for nBody=1:NmechBodies
    Cval{nBody}= zeros(nSparse,1) ;
    Cvalpre{nBody}= zeros(nSparse,1);
end
CvalReg         = zeros(nSparse,1);
CvalpreReg         = zeros(nSparse,1);
Mval            = zeros(nSparse,1);
Kdval              = [];
Cdval              = [];
CdvalReg              = [];
Mdval              = [];
Idir             = [];
Jdir              = [];
% Initialise variables
count           = 0;
nmst            = 0;
nmstDir         = 0;

FlagMechVolume=Unknown.FlagMechVolume;
Ccond=cell(NmechBodies,1);
Ccondpre=cell(NmechBodies,1);
%=========================================================================
% Loop through elements
%=========================================================================
for i=1:nelem
    for nBody=1:NmechBodies
    Ccond{nBody}=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    Ccondpre{nBody}=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    end
    CReg=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    CpreReg=zeros(esizet+3*esizeH1,esizet+3*esizeH1);
    kappalowV=zeros(3,1);
    flagCond=zeros(3,1);
    
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
    u(lunkvM>0)=X(lunkvM(lunkvM>0));
    
    % Store the dynamic solution in the array
    Xdir  = [A;u];
    
    %---------------------------------------------------------------------
    % Determine the Elemental matrices
    %---------------------------------------------------------------------
    [M,C,K,Res]=GaussInteg(Basis,Quadrature,ProblemData,...
        flag,mue,Adc,A0,mycoord,eltype,probstatic,material,subFlag(i),...
        i,lec,lfc,xy,cond,coupling,FlagVolume,SourceMapping);
    
    
    %-------------------------------------------------------------------------
    % Define kappa coefficient defined by regularization term and conductivity
    %-------------------------------------------------------------------------
    
    if subFlag(i)==2                                   % Conducting region
        if probstatic==0                            % Dynamic problem
            kappalow=sigma(mat(i))*muz;
            kappalowReg=0;
            for nBody=1:NmechBodies
                if any(ismember(mat(i),matMech(nBody,:)))
                    kappalowV(nBody)=kappalow;
                    flagCond(nBody)=kappalow;
                end
            end
            
        else                       % Static problem
            kappalowV=zeros(3,1);
            kappalowReg=regTerm;
            flagCond=zeros(3,1);
        end
    else
        flagCond=zeros(3,1);
        kappalowV=zeros(3,1);
        kappalowReg=regTerm;
    end
    
    % Redefine regularisation for solver option 1
    if ProblemData.sol.regopt ==1
        if subFlag(i)~=2
           kappaV=zeros(3,1);
           kappaReg=0;
        else
            kappaV=kappalowV;
            kappaReg=kappalowReg;
        end
    else
        kappaV=kappalowV;
        kappaReg=kappalowReg;
    end
    
    
    % Regularization
    for nBody=1:NmechBodies
    Ccond{nBody}(1:6,1:6)=kappalowV(nBody)*C(1:6,1:6);
    Ccondpre{nBody}(1:6,1:6)=Ccond{nBody}(1:6,1:6);
    end
    CReg(1:6,1:6)=kappalowReg*C(1:6,1:6);
    CpreReg(1:6,1:6)=CReg(1:6,1:6);
    
    % Extract number of interior EM basis
    nintbas=Unknown.EM.nintbas;
    
    
    %  Add regularisation to higher order block
    for p=1:esizet
        if p<= 6
            for nBody=1:NmechBodies
            Ccond{nBody}(p,7:esizet)=kappaV(nBody)*C(p,7:esizet);
            end
            CReg(p,7:esizet)=kappaReg*C(p,7:esizet);
        else
            for nBody=1:NmechBodies
            Ccond{nBody}(p,1:esizet)=kappaV(nBody)*C(p,1:esizet);
            Ccondpre{nBody}(p,1:esizet)=abs(kappaV(nBody))*C(p,1:esizet);
            end
            CReg(p,1:esizet)=kappaReg*C(p,1:esizet);
            CpreReg(p,1:esizet)=abs(kappaReg)*C(p,1:esizet);
        end
    end
    
    if ProblemData.sol.regopt==1
        % high order edges (gradients)
        for nBody=1:NmechBodies
        Ccondpre{nBody}(7:6*(order+1),7:6*(order+1))=abs(kappaV(nBody))*C(7:6*(order+1),7:6*(order+1));
        end
        CpreReg(7:6*(order+1),7:6*(order+1))=abs(kappaReg)*C(7:6*(order+1),7:6*(order+1));
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappaV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappaReg)*C(st:en,st:en);
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappalowV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappaV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappaReg)*C(st:en,st:en);
        
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappalowV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
    else
        % high order edges (gradients)
        for nBody=1:NmechBodies
        Ccondpre{nBody}(7:6*(order+1),7:6*(order+1))=abs(kappalowV(nBody))*C(7:6*(order+1),7:6*(order+1));
        end
        CpreReg(7:6*(order+1),7:6*(order+1))=abs(kappalowReg)*C(7:6*(order+1),7:6*(order+1));
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappalowV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappalowV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappalowV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        for nBody=1:NmechBodies
        Ccondpre{nBody}(st:en,st:en)=abs(kappalowV(nBody))*C(st:en,st:en);
        end
        CpreReg(st:en,st:en)=abs(kappalowReg)*C(st:en,st:en);
    end
    
    
    % Complete damping matrix (no regularisation for mechanic block)
    for nBody=1:NmechBodies
    Ccond{nBody}(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1)=flagCond(nBody)*C(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1);
    Ccondpre{nBody}(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1)=Ccond{nBody}(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1);
    Ccond{nBody}(1:esizet,esizet+1:esizet+3*esizeH1)=flagCond(nBody)*C(1:esizet,esizet+1:esizet+3*esizeH1);
    Ccondpre{nBody}(1:esizet,esizet+1:esizet+3*esizeH1)=Ccond{nBody}(1:esizet,esizet+1:esizet+3*esizeH1);
    end
    %---------------------------------------------------------------------
    % Linear system assembly in vector format
    %---------------------------------------------------------------------
    [Mval,Cval,Cvalpre,CvalReg,CvalpreReg,Kval,Resid,I,J,Idir,Jdir,nmst,nmstDir,Kdval,Cdval,CdvalReg,Mdval]=matrixAssembly(K,Ccond,Ccondpre,CReg,CpreReg,M,Res,Xdir,Mval,Cval,Cvalpre,CvalReg,CvalpreReg,Kval,Resid,I,J,Idir,Jdir,lunkv,probstatic,esizet,esizeH1,nmst,nmstDir,omega,Kdval,Cdval,CdvalReg,Mdval,freqSweep,NmechBodies);
    
end


%=========================================================================
% Construct the Linear system in sparse format
%=========================================================================
% Build stiffness and damping matrices from I,J,Kval and Cval
stiff  = sparse(I(1:nmst),J(1:nmst),Kval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damp=cell(NmechBodies,1);
damppre=cell(NmechBodies,1);
for nBody=1:NmechBodies
damp{nBody}  = sparse(I(1:nmst),J(1:nmst),Cval{nBody}(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damppre{nBody}=sparse(I(1:nmst),J(1:nmst),Cvalpre{nBody}(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
end
dampReg  = sparse(I(1:nmst),J(1:nmst),CvalReg(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
mass  = sparse(I(1:nmst),J(1:nmst),Mval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damppreReg=sparse(I(1:nmst),J(1:nmst),CvalpreReg(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
Kd=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),Kdval(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
% Cd=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),Cdval(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
Cd=[];
CdReg=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),CdvalReg(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
Md=sparse(Idir(1:nmstDir),Jdir(1:nmstDir),Mdval(1:nmstDir),nunktEM+nunktMech,ndirEM+ndirM);
end


