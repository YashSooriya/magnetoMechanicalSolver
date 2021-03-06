function probdata = problemMHIGradXNewCoilOnlyEnds(pm)


% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Choose the BC type for the outer boundary
probdata.bctypeouter=2;   % 2 Dirichlet outer BC
                          % 3 Neumann outer BC

% Polynomial Degree Info Order
order=pm;

%--------------------------------------------------------------------------
% Job data
%--------------------------------------------------------------------------
% job=sprintf('TestGradXUnionBackground_SplitPartitions');
job=sprintf('MHIGradXNewCoil2FixedOnlyEnds');
% job='sphere1';3
meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
                   % 4= ansys, 5 = Opera

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;

%--------------------------------------------------------------------------
% Define solver option
%--------------------------------------------------------------------------


probdata.sol.regopt =2;     % 1- use regularisation with cg for gradent blocks
                            % 2- use regularisation with direct solve for gradient blocks
                            % 3- direct solve
                            
% Define tolerance for GMRES solver
TOL_GMRES=1e-5;
probdata.TOL_GMRES=TOL_GMRES;

% -------------------------------------------------------------------------
% Define current density magnitude (DC and AC coils)
JDC         = 253.2*1e6;
Jratio      = 0.022116904;
JAC         = JDC*Jratio;

%--------------------------------------------------------------------------
% Material properties
%--------------------------------------------------------------------------

nmat = 10;                            % Number of materials
muz = 1.256637061435917e-06;          % Mu_z (Free sopace permeability)
epz = 0;                              % Ep_z
jsrc=[0 0 0 0 0 0];
% Conductivity, gamma
regTerm=1e-4;
sigma   = [71e6 71e6 71e6 33e6 33e6 33e6 1.4e6 1.4e6 1.4e6 regTerm regTerm regTerm];

% Density, rho
rho     = [2710 2710 2710 2698 2698 2698 7900 7900 7900 0 0 0];

% Young's Modulus, E
E       = [81 81 81 81 81 81 210 210 210 0 0 0]*1e9;
% Mu_r, Ep_r, Sigma, J

% Poisson's Ratio, nu
nu      = [0.338 0.338 0.338 0.337 0.337 0.337 0.283 0.283 0.283 0 0 0];

mu_r=[1 1 1 1 1 1 1 1 1 1 1 1];
mu=mu_r;

% Define Lame constants
lambda   = (E.*nu)./((1+nu).*(1-2.*nu));
G        = E./(2.*(1+nu));

% For this problem we have two
D_1=zeros(6);
D_1(1,1)=lambda(1)+2*G(1);
D_1(2,2)=lambda(1)+2*G(1);
D_1(3,3)=lambda(1)+2*G(1);
D_1(4,4)=G(1);
D_1(5,5)=G(1);
D_1(6,6)=G(1);
D_1(1,2)=lambda(1);
D_1(1,3)=lambda(1);
D_1(2,1)=D_1(1,2);
D_1(3,1)=D_1(1,3);
D_1(2,3)=lambda(1);
D_1(3,2)=D_1(2,3);


D_2=zeros(6);
D_2(1,1)=lambda(4)+2*G(4);
D_2(2,2)=lambda(4)+2*G(4);
D_2(3,3)=lambda(4)+2*G(4);
D_2(4,4)=G(4);
D_2(5,5)=G(4);
D_2(6,6)=G(4);
D_2(1,2)=lambda(4);
D_2(1,3)=lambda(4);
D_2(2,1)=D_2(1,2);
D_2(3,1)=D_2(1,3);
D_2(2,3)=lambda(4);
D_2(3,2)=D_2(2,3);


D_3=zeros(6);
D_3(1,1)=lambda(7)+2*G(7);
D_3(2,2)=lambda(7)+2*G(7);
D_3(3,3)=lambda(7)+2*G(7);
D_3(4,4)=G(7);
D_3(5,5)=G(7);
D_3(6,6)=G(7);
D_3(1,2)=lambda(7);
D_3(1,3)=lambda(7);
D_3(2,1)=D_3(1,2);
D_3(3,1)=D_3(1,3);
D_3(2,3)=lambda(7);
D_3(3,2)=D_3(2,3);

D_4=zeros(6);

% Store the tensor in a cell array
D{1}=D_1;
D{2}=D_1;
D{3}=D_1;
D{4}=D_2;
D{5}=D_2;
D{6}=D_2;
D{7}=D_3;
D{8}=D_3;
D{9}=D_3;
D{10}=D_4;
D{11}=D_4;
D{12}=D_4;
% Specify the material to be used a conductor (where gradients basis
% functions to be included)
% Also specify mechanical subdomain

matcond=[1 2 3 4 5 6 7 8 9 12];
subdom_mech=[1 2 3 4 5 6 7 8 9];
mur=mu(2)/mu(1);

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 1e-3;  % Object size (reescale due to units issue)
shift=[0 0 0]; % Object shift

%-------------------------------------------------------------------------------
% Define material numbers of conducting shields to be used in DOF numbering
%-------------------------------------------------------------------------------
mat4K=[1 2 3];
mat77K=[4 5 6];
matOVC=[7 8 9];
matAir=10;
%-------------------------------------------------------------------------------
% Store material properties on probdata.matr substructure
%-------------------------------------------------------------------------------

probdata.matr.muz=muz;
probdata.matr.mur=mur;
probdata.matr.epz=epz;
probdata.matr.mu=mu;
probdata.matr.sigma=sigma;
probdata.matr.jsrc=jsrc;
probdata.matr.delta=delta;
probdata.matr.shift=shift;
probdata.matr.matcond=matcond;
probdata.matr.subdom_mech=subdom_mech;
probdata.matr.regTerm=regTerm;
probdata.matr.rho=rho;
probdata.matr.E=E;
probdata.matr.nu=nu;
probdata.matr.lambda=lambda;
probdata.matr.G=G;
probdata.matr.D=D;
probdata.mat4K=mat4K;
probdata.mat77K=mat77K;
probdata.matOVC=matOVC;
probdata.matAir=matAir;

%-------------------------------------------------------------------------
% Blending Function Info (over write defaults)
%-------------------------------------------------------------------------
gorder =1;        % order (if gorder > 0, quadlin =2 required)
g1 = 2;           % use exact geometry for a sphere

sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once
svchk = 0;         % Check surface/volume as expecting a sphere


probdata.jb.gorder = gorder;
probdata.jb.g1 = g1;
probdata.jb.sufv = sufv(1);
probdata.jb.gag = gag(1);
probdata.mesh.svchk=svchk;

% Define input argument for functions
arg.sigma=sigma;
arg.mu=mu(2)*muz;
arg.muz=muz;
arg.mur=mur;
arg.alpha=delta;
arg.muz=muz;
arg.nu=nu;
arg.E=E;
arg.lambda=lambda;
arg.G=G;
arg.JDC=JDC;
arg.JAC=JAC;


% define the function handle for Dirichlet BC's
probdata.es.dirfun=@esproblemdir;
probdata.es.dirfunarg=arg;

% define the source terms
probdata.es.srcfunarg=[];
probdata.es.srcfun=@esproblemsrc;
probdata.es.srcfunarg=arg;

% define the source terms
probdata.es.srcfunargmech=[];
probdata.es.srcfunmech=@esproblemsrcmech;
probdata.es.srcfunargmech=arg;

% define the function handle for Neumann BC's
probdata.es.neufun=@esproblemneu;
probdata.es.neufunarg=arg;

% exact solution
probdata.es.exactfun=@esproblemexact;
probdata.es.exactfunarg=arg;


% Gradient of exact solution
probdata.es.GradientExact=@gradExact;
probdata.es.gradfunarg=arg;

% exact magnetic field
probdata.es.exactcurlfun=@esproblemexactcurl;
probdata.es.exactcurlfunarg=arg;



% Dirichlet Boundary Condition----------------------------------------------
function e=esproblemdir(x,y,z,index,arg,probstatic,probFlag,omega)

if probFlag==1                  % EM problem
    e=[0;0;0];
    
    
elseif probFlag==2                      % Mechanical problem
    
    e=[0;0;0];
end

function grade=gradExact(x,y,z,arg)

% define the gradient components
grade=zeros(3,3);
grade(1,1)=0;
grade(1,2)=0;
grade(1,3)=0;
grade(2,1)=0;
grade(2,2)=0;
grade(2,3)=0;
grade(3,1)=0;
grade(3,2)=0;
grade(3,3)=0;


% end

% Neumann Boundary Condition-----------------------------------------------

function curle=esproblemneu(x,y,z,index,arg,probstatic)
curle=zeros(3,1);

% Src term-----------------------------------------------------------------

function src=esproblemsrc(x,y,z,arg,domain,probstatic)               % Source term for es
phi=atan2(y,x);
JDC=arg.JDC;
JAC=arg.JAC;
R=0.1979;
y1=-0.13335;
y2=0.13335;
x2=sqrt(R^2-y2^2);

if domain==11
    if probstatic==1
        src=1.256637061435917e-06*[-sin(phi)*JDC; cos(phi)*JDC; 0];
    else
        src=[0;0;0];
    end
    


elseif domain==12
    z1=0.1629;


    if probstatic==0 || probstatic==2
        
        
        if y>y1 && y<y2
            if z<z1
                src=-1.256637061435917e-06*JAC*[-sin(phi); cos(phi); 0];
            elseif z>z1
                src=1.256637061435917e-06*JAC*[-sin(phi); cos(phi); 0];
            end
        elseif y>=y2 
        zp=z-z1;
         phic=atan2(y2,x2);
        yp=R*(phi-phic);
        theta=atan2(zp,yp);

           dirVec=-[sin(phi)*sin(theta); -cos(phi)*sin(theta); cos(theta)];
            src=1.256637061435917e-06*JAC*dirVec/norm(dirVec);

        elseif y<=y1
        zp=z-z1;
         phic=atan2(y1,x2);
        yp=R*(phi-phic);
        theta=atan2(zp,yp);

            dirVec=-[sin(phi)*sin(theta); -cos(phi)*sin(theta); cos(theta)];
            src=1.256637061435917e-06*JAC*dirVec/norm(dirVec);
        end
    else
        src=[0; 0; 0];
    end
    
else
    src=[0; 0; 0];
end

function srcmech=esproblemsrcmech(x,y,z,dum,arg,domain) % Source term for es
srcmech=[0;0;0];

% Analytical Solution Field------------------------------------------------
function [out]=esproblemexact(x,y,z,arg,probstatic,probFlag,omega)

out=zeros(3,1);


% Analytical Solution for B=curl A =mu H-----------------------------------
% Gives total field as output
function [curle,hf,H0]=esproblemexactcurl(x,y,z,arg,probstatic,omega)
curle=zeros(3,1);
hf=zeros(3,1);
H0=zeros(3,1);




