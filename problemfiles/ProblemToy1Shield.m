function probdata = ProblemToy1Shield(pm,option)


% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Choose the BC type for the outer boundary
probdata.bctypeouter=2;   % 2 Dirichlet outer BC
                          % 3 Neumann outer BC

% Polynomial Degree Info Order
order=pm; 
opt=option;

%--------------------------------------------------------------------------
% Define postprocessing switchs
%--------------------------------------------------------------------------

% % plot the perturbed solution on a line
% % plotoption = 0 - do not plot
% % plotoption = 1 - plot
 plotoption = 1;
 
% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 0;
 
% Calculate the H(curl) norm of the error
% % erroroption = 0 - do not calculate error
% % erroroption = 1 - calculate error
 erroroption = 1; % Exact E field not entered!
 
% % output the VTK file
% % vtkoption = 0 - do not output
% % vtkoption = 1 - output
 vtkoption =0;
 
% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 0; 

% % plot the  H_alpha-H_0 - D2G PT H_0  on a line
% % darkooption = 0 - do not plot
% % darkooption = 1 - plot
 darkooption = 0; 
 
probdata.mesh.plotoption=plotoption;
probdata.mesh.erroroption=erroroption;
probdata.mesh.vtkoption=vtkoption;
probdata.mesh.voltoption=voltoption;
probdata.mesh.darkooption=darkooption;

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

%job=sprintf('sphere%d_2',option);   % Job Name 
%job=sprintf('Toy_noShields');
%job=sprintf('Toy_noShields2');
%job=sprintf('Toy_noGrad3');
job=sprintf('Toy_1shield');


% job='sphere1';3
meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
                   % 4= ansys, 5 = Opera
%probdata.jb.name=sprintf('sphere_small_%d',option);
probdata.jb.name=sprintf('Opera_sphere_mesh_lin_100');
probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;

%--------------------------------------------------------------------------
% Define solver option
%--------------------------------------------------------------------------


probdata.sol.regopt =3;     % 1- use regularisation with cg for gradent blocks
                            % 2- use regularisation with direct solve for gradient blocks
                            % 3- direct solve
                            
% -------------------------------------------------------------------------
% Define current magnitude (DC and AC coils)
% Define the static current density magnitude
JDC         = 253.2*1e6;
Jratio=0.022116904;
JAC         = JDC*Jratio;
                            
 %--------------------------------------------------------------------------
% Material properties
%--------------------------------------------------------------------------

nmat = 4;                            % Number of materials
muz = 1.256637061435917e-06;         % Mu_z (Free sopace permeability)
epz = 0;                             % Ep_z
%omega = 2*3.141592*50;               % Omega (angular frequency)
omega=133.5;
jsrc=[0 0 0 0];
% Conductivity, gamma
sigma   = [0 0 71 71 0]*1e6;
regTerm=1e-1;

% Density, rho
rho     = [0 0 2710 2710 0];

% Young's Modulus, E
E       = [0 0 81 81 0]*1e9;
% Mu_r, Ep_r, Sigma, J

% Poisson's Ratio, nu
nu      = [0 0 0.338 0.338 0];

mu_r=[1 1 1 1 1];
mu=mu_r;

% Define Lame constants
lambda   = (E.*nu)./((1+nu).*(1-2.*nu));
G        = E./(2.*(1+nu)); 

% For this problem we have two
D_1=zeros(6);
D_1(1,1)=lambda(3)+2*G(3);
D_1(2,2)=lambda(3)+2*G(3);
D_1(3,3)=lambda(3)+2*G(3);
D_1(4,4)=G(3);
D_1(5,5)=G(3);
D_1(6,6)=G(3);
D_1(1,2)=lambda(3);
D_1(1,3)=lambda(3);
D_1(2,1)=D_1(1,2);
D_1(3,1)=D_1(1,3);
D_1(2,3)=lambda(3);
D_1(3,2)=D_1(2,3);

D_2=zeros(6);



% Store the tensor in a cell array
D{1}=D_2;
D{2}=D_2;
D{3}=D_1;
D{4}=D_1;
D{5}=D_2;



% Specify the material to be used a conductor (where gradients basis
% functions to be included)
% Also specify mechanical subdomain

matcond=[1 2 3 4];
subdom_mech=[3 4];
mur=mu(2)/mu(1);

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 1;  % Object size
shift=[0 0 0]; % Object shift

%-------------------------------------------------------------------------------
% Store material properties on probdata.matr substructure
%-------------------------------------------------------------------------------

probdata.matr.muz=muz;
probdata.matr.mur=mur;
probdata.matr.epz=epz;
probdata.matr.omega=omega;
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

%-------------------------------------------------------------------------
% Blending Function Info (over write defaults)
%-------------------------------------------------------------------------
gorder =2;        % order (if gorder > 0, quadlin =2 required)
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
p_0=1e9;
arg.sigma=sigma(2);
arg.omega=omega;
arg.mu=mu(2)*muz;
arg.muz=muz;
arg.mur=mur;
arg.alpha=delta;
arg.muz=muz;
arg.nu=nu;
arg.E=E;
arg.p_a=p_0;
arg.p_b=p_0*1000;
arg.lambda=lambda(3);
arg.G=G(3);
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
function e=esproblemdir(x,y,z,index,arg,probstatic,probFlag)    

if probFlag==1                  % EM problem
e=[0;0;0];


elseif probFlag==2                      % Mechanical problem
    
    e=[x^2+y^2+z^2;0;0];
end
 
    function grade=gradExact(x,y,z,arg)

% define the gradient components
grade=zeros(3,3);
grade(1,1)=2*x;
grade(1,2)=2*y;
grade(1,3)=2*z;
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

function src=esproblemsrc(x,y,z,dum,arg,domain,probstatic)               % Source term for es
    phi=atan2(y,x);
    JDC=arg.JDC;
    JAC=arg.JAC;
if domain==1 || domain==2
    src=1.25e-6*[-sin(phi)*JDC; cos(phi)*JDC; 0];
 elseif domain==3 || domain==4
     if probstatic==0
    src=[-sin(phi)*JAC; cos(phi)*JAC; 0];
     else
         src=[0;0;0];
     end
else 
src=[0;0;0];
end

function srcmech=esproblemsrcmech(x,y,z,dum,arg,domain) % Source term for es
    lambda=arg.lambda;
    G=arg.G;
srcmech=[2*lambda+8*G;0;0];

% Analytical Solution Field------------------------------------------------
function [out]=esproblemexact(x,y,z,arg,probstatic,probFlag)   
    if probFlag==1

out=zeros(3,1);

    elseif probFlag==2
        out=[x^2+y^2+z^2;0;0];
    end

 
% Analytical Solution for B=curl A =mu H-----------------------------------
% Gives total field as output
function [curle,hf,H0]=esproblemexactcurl(x,y,z,arg,probstatic)  
curle=zeros(3,1);
hf=zeros(3,1);
H0=zeros(3,1);




