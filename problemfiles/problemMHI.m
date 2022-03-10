function probdata = problemMHI(pm)


% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Choose the BC type for the outer boundary
probdata.bctypeouter=2;   % 2 Dirichlet outer BC
                          % 3 Neumann outer BC

% Polynomial Degree Info Order
order=pm; 


 %--------------------------------------------------------------------------
 % Define lines for plotting the magnetic fields
 %--------------------------------------------------------------------------
 
% %interpolation number
 N = 200;
 probdata.mesh.N = N;
 N1=200;


 starp = [0.0,0,-1.19];   % startpoin of the line
 overp =[0.0,0,1.19];   % endpoint of the line
      point = starp;   
     deltax = (overp(1)-starp(1))/(N1-1);
     deltay = (overp(2)-starp(2))/(N1-1);
     deltaz = (overp(3)-starp(3))/(N1-1);

     for i = 2:N1
         point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
     end 

     
 probdata.mesh.point = point;

%--------------------------------------------------------------------------
% Job data
%--------------------------------------------------------------------------

job=sprintf('MHI');       % Mesh name

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
% Define current magnitude (DC and AC coils)
% Define the static current density magnitude
JDC         = 253.2*1e6;
Jratio      = 0.022116904;
JAC         = JDC*Jratio;
                            
 %--------------------------------------------------------------------------
% Material properties
%--------------------------------------------------------------------------

nmat = 8;                            % Number of materials
muz = 1.256637061435917e-06;         % Mu_z (Free sopace permeability)
epz = 0;                             % Ep_z
jsrc=[0 0 0 0 0 0];
% Conductivity, gamma
regTerm=1e-4;
sigma   = [regTerm regTerm regTerm regTerm 71e6 33e6 1.4e6 regTerm];

% Density, rho
rho     = [0 0 0 0 2710 2698 7900 0];

% Young's Modulus, E
E       = [0 0 0 0 81 81 210 0]*1e9;
% Mu_r, Ep_r, Sigma, J

% Poisson's Ratio, nu
nu      = [0 0 0 0 0.338 0.337 0.283 0];

mu_r=[1 1 1 1 1 1 1 1];
mu=mu_r;

% Define Lame constants
lambda   = (E.*nu)./((1+nu).*(1-2.*nu));
G        = E./(2.*(1+nu)); 

% For this problem we have two
D_1=zeros(6);
D_1(1,1)=lambda(5)+2*G(5);
D_1(2,2)=lambda(5)+2*G(5);
D_1(3,3)=lambda(5)+2*G(5);
D_1(4,4)=G(5);
D_1(5,5)=G(5);
D_1(6,6)=G(5);
D_1(1,2)=lambda(5);
D_1(1,3)=lambda(5);
D_1(2,1)=D_1(1,2);
D_1(3,1)=D_1(1,3);
D_1(2,3)=lambda(5);
D_1(3,2)=D_1(2,3);


D_2=zeros(6);
D_2(1,1)=lambda(6)+2*G(6);
D_2(2,2)=lambda(6)+2*G(6);
D_2(3,3)=lambda(6)+2*G(6);
D_2(4,4)=G(6);
D_2(5,5)=G(6);
D_2(6,6)=G(6);
D_2(1,2)=lambda(6);
D_2(1,3)=lambda(6);
D_2(2,1)=D_2(1,2);
D_2(3,1)=D_2(1,3);
D_2(2,3)=lambda(6);
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
D{1}=D_4;
D{2}=D_4;
D{3}=D_4;
D{4}=D_4;
D{5}=D_1;
D{6}=D_2;
D{7}=D_3;
D{8}=D_1;

% Specify the material to be used a conductor (where gradients basis
% functions to be included)
% Also specify mechanical subdomain

matcond=[3 4 5 6 7];
subdom_mech=[5 6 7];
mur=mu(2)/mu(1);

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 1;     % Object size
shift=[0 0 0]; % Object shift

%-------------------------------------------------------------------------------
% Define material numbers of conducting shields to be used in DOF numbering
%-------------------------------------------------------------------------------
mat4K=[5];
mat77K=[6];
matOVC=[7];
matAir=8;

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
    if z<0
        JAC=-JAC;
    end
if domain==1 || domain==2
    if probstatic==1
    src=1.256637061435917e-06*[-sin(phi)*JDC; cos(phi)*JDC; 0];
    else
        src=[0;0;0];
    end
 elseif domain==3 || domain==4
     if probstatic==0
    src=1.256637061435917e-06*[-sin(phi)*JAC; cos(phi)*JAC; 0];
     else
         src=[0;0;0];
     end
else 
src=[0;0;0];
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



