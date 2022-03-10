function [M,C,K,Res,Volume,AreaInt]=GaussInteg2(Basis,DirichletFaceBasis,Quadrature,probdata,flag,gorder,mue,esizet,esizeH1,omega,Adc,Xa,Xu,mycoord,eltype,probstatic,material,subFlag,ielem,lec,lfc,xy,cond,matc,Volume,AreaInt)
%=========================================================================
% Determine the elemental stiffness matrix for triangles
%=========================================================================
% Extract the relevant data from the input structures
if eltype==1
    ephx=Basis.ephx1;
    ephy=Basis.ephy1;
    ephz=Basis.ephz1;
    ecphx=Basis.ecphx1;
    ecphy=Basis.ecphy1;
    ecphz=Basis.ecphz1;
    gphx=Basis.gphx1;
    gphy=Basis.gphy1;
    gphz=Basis.gphz1;
    phh1=Basis.phh11;
    H1basis=Basis.H1basis1;
    gradH1basisx=Basis.gradH1basis1x;
    gradH1basisy=Basis.gradH1basis1y;
    gradH1basisz=Basis.gradH1basis1z;
    ephdfx=DirichletFaceBasis.ephdfx1;
ephdfy=DirichletFaceBasis.ephdfy1;
ephdfz=DirichletFaceBasis.ephdfz1;
ecphdfx=DirichletFaceBasis.ecphdfx1;
ecphdfy=DirichletFaceBasis.ecphdfy1;
ecphdfz=DirichletFaceBasis.ecphdfz1;
phh1f=DirichletFaceBasis.phh1f1;
phH1f=DirichletFaceBasis.phH1f1;
gphfx=DirichletFaceBasis.gphfx1;
gphfy=DirichletFaceBasis.gphfy1;
gphfz=DirichletFaceBasis.gphfz1;
elseif eltype==2
      ephx=Basis.ephx2;
    ephy=Basis.ephy2;
    ephz=Basis.ephz2;
    ecphx=Basis.ecphx2;
    ecphy=Basis.ecphy2;
    ecphz=Basis.ecphz2;
    gphx=Basis.gphx2;
    gphy=Basis.gphy2;
    gphz=Basis.gphz2;
    phh1=Basis.phh12;
    H1basis=Basis.H1basis2;
    gradH1basisx=Basis.gradH1basis2x;
    gradH1basisy=Basis.gradH1basis2y;
    gradH1basisz=Basis.gradH1basis2z;
        ephdfx=DirichletFaceBasis.ephdfx2;
ephdfy=DirichletFaceBasis.ephdfy2;
ephdfz=DirichletFaceBasis.ephdfz2;
ecphdfx=DirichletFaceBasis.ecphdfx2;
ecphdfy=DirichletFaceBasis.ecphdfy2;
ecphdfz=DirichletFaceBasis.ecphdfz2;
phh1f=DirichletFaceBasis.phh1f2;
phH1f=DirichletFaceBasis.phH1f2;
gphfx=DirichletFaceBasis.gphfx2;
gphfy=DirichletFaceBasis.gphfy2;
gphfz=DirichletFaceBasis.gphfz2;
end

intfw=Quadrature.intfw;
intfxi=Quadrature.intfxi;
intfet=Quadrature.intfet;
nipf=Quadrature.nipf;


%=========================================================================
% Initialisation of the linear system vectors/matrices
%=========================================================================

% Predefine the size of the blocks
R_A        = zeros(esizet,1);
R_U        = zeros(3*esizeH1,1);
R_UB        = zeros(3*esizeH1,1);

% Predefine the size of the blocks
K_AA       = zeros(esizet);
C_AA       = zeros(esizet);
K_AU       = zeros(esizet,3*esizeH1);
K_UA       = zeros(3*esizeH1,esizet);
C_AU       = zeros(esizet,3*esizeH1);
C_UA       = zeros(3*esizeH1,esizet);
M_AA       = zeros(esizet);
M_UU       = zeros(3*esizeH1);
M_AU       = zeros(esizet,3*esizeH1);
M_UA       = zeros(3*esizeH1,esizet);
C_UU       = zeros(3*esizeH1);



% Extract the elemental values
A     = Xa(1:esizet,1);
u     = Xu(1:3*esizeH1,1);




% Define the DC field as the magnetic field for non linearity
%if (probstatic == 1)
    %Adc = A;
%end

%--------------------------------------------------------------------------
% Define the elasticity tensor C_ijkl=D
%--------------------------------------------------------------------------

if subFlag==2
D=probdata.matr.D{material};
else
    D=0;
end
lambda=probdata.matr.lambda(material);
G=probdata.matr.G(material);

%=========================================================================
% Compute the body terms
%=========================================================================

    
    
    %---------------------------------------------------------------------
    % Compute the body forces
    %---------------------------------------------------------------------
    % Determine the current source value

 

    
    %-----------------------------------------------------------------
    % Construct the Directional derivative blocks
    %-----------------------------------------------------------------
    % Compute the linearised mass matrix components
    [M_UU]=massLinear(Quadrature,gphx,gphy,gphz,H1basis,esizeH1,flag,gorder,mycoord,subFlag,probdata,material);
    
    % Compute the linearised damping matrix components
    [C_AA,C_UU]=dampLinear(Quadrature,ephx,ephy,ephz,gphx,gphy,gphz,esizet,esizeH1,flag,gorder,mycoord);
    
    % Compute the linearised stiffness matrix components
    [K_AA,K_AU,K_UA,K_UU]=stiffLinear(Quadrature,ecphx,ecphy,ecphz,gphx,gphy,gphz,gradH1basisx,gradH1basisy,gradH1basisz,esizet,esizeH1,flag,gorder,mue,mycoord,D,subFlag,A);
    %[K_AA,K_UU]=stiffLinear1(Quadrature,ecphx,ecphy,ecphz,gphx,gphy,gphz,gradH1basisx,gradH1basisy,gradH1basisz,esizet,esizeH1,flag,gorder,mue,mycoord,lambda,G,subFlag);
    [K_UA]=stiffLinearBoun(K_UA,xy,esizet,esizeH1,intfxi,intfet,intfw,nipf,eltype,cond,...
    ielem,lec,lfc,flag,gorder,probdata,...
    mycoord,phh1f,phH1f,gphfx,gphfy,gphfz,ecphdfx,ecphdfy,ecphdfz,mue,probstatic,A,matc);
    
    % Construct the Residual vector components
    [R_A,R_U,Volume]=residual(R_A,R_U,Quadrature,ephx,ephy,ephz,ecphx,ecphy,ecphz,gphx,gphy,gphz,phh1,gradH1basisx,gradH1basisy,gradH1basisz,H1basis,esizet,esizeH1,flag,gorder,mue,mycoord,A,u,probdata,material,probstatic,D,subFlag,Volume);
    [R_U,AreaInt]= residualBoun2(R_U,xy,esizet,esizeH1,intfxi,intfet,intfw,nipf,eltype,cond,...
    ielem,lec,lfc,flag,gorder,probdata,...
    mycoord,phh1f,phH1f,gphfx,gphfy,gphfz,ecphdfx,ecphdfy,ecphdfz,mue,probstatic,A,matc,AreaInt);


%=========================================================================
% Compute the edge/boundary terms
%=========================================================================

% Not necessary for purely electromganetics (construct for mechanics)

%=========================================================================
% Construct the local (element-wise) Linear system 
%=========================================================================
% Add Rayleigh damping term to the linearised mechanical damping matrix


% Construct the linearised mass matrix
M    = [M_AA M_AU;M_UA M_UU];

% Construct the linearised damping matrix
C    = [C_AA C_AU;C_UA C_UU];

% Construct the linearised stiffness matrix
K    = [K_AA K_AU;K_UA K_UU];

% Local residual vector
Res  = [R_A;R_U];

end