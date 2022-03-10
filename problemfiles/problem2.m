function probdata=problem2(order,option,Jrat)

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order

% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Choose the BC type for the outer boundary
probdata.bctypeouter=3;   % 2 Dirichlet outer BC
                          % 3 Neumann outer BC


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
 N = 100;
 probdata.mesh.N = N;
 N1=10;
 N2=N-N1;

 starp = [0.0,0,0.00];   % startpoin of the line
 overp =[0.0,0,0.009];   % endpoint of the line
      point = starp;   
     deltax = (overp(1)-starp(1))/(N1-1);
     deltay = (overp(2)-starp(2))/(N1-1);
     deltaz = (overp(3)-starp(3))/(N1-1);

     for i = 2:N1
         point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
     end 
 starp = [0.0,0,0.011];   % startpoin of the line
 overp =[0.0,0,1];   % endpoint of the line
     
      point = [point; starp];   
     deltax = (overp(1)-starp(1))/(N2-1);
     deltay = (overp(2)-starp(2))/(N2-1);
     deltaz = (overp(3)-starp(3))/(N2-1);

     for i = 2:N2
         point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
     end 
     
     
 probdata.mesh.point = point;


%--------------------------------------------------------------------------
% Job data
%--------------------------------------------------------------------------

job=sprintf('sphere%d_2',option);   % Job Name 



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


probdata.sol.regopt =2;     % 1- use regularisation with cg for gradent blocks
                            % 2- use regularisation with direct solve for gradient blocks
                            % 3- direct solve
                            
%--------------------------------------------------------------------------
% Material properties
%--------------------------------------------------------------------------

nmat = 2;                            % Number of materials
muz = 1.256637061435917e-06;         % Mu_z (Free sopace permeability)
epz = 0;                             % Ep_z
freq=50;
omega = 2*3.141592*freq;               % Omega (angular frequency)
% Mu_r, Ep_r, Sigma, J

% Mat 1
mu(1) = 1;                           % Permeability (relative)
epl(1) = 0;                          % Permittivity (relative)
sigma(1) = 1e-1;                     % Conductivity
regTerm=1e-1;                        % Regularization term
% Mechanical properties
rho(1) = 7800;                       % Density
E(1) = 210e9;                        % Young's modulus
nu(1) = 0.33;                        % Poisson's ratio

% Mat 2
mu(2) = 1;                         % Permeability (relative)
epl(2) = 0;                          % Permittivity (relative)
sigma(2) = 6e7;                   % Conductivity
% Mechanical properties
rho(2)=0;                            % Density
E(2)=0;                              % Young's modulus
nu(2)=0;                             % Poisson's ratio

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

D_2(1,1)=lambda(2)+2*G(2);
D_2(2,2)=lambda(2)+2*G(2);
D_2(3,3)=lambda(2)+2*G(2);
D_2(4,4)=G(2);
D_2(5,5)=G(2);
D_2(6,6)=G(2);
D_2(1,2)=lambda(2);
D_2(1,3)=lambda(2);
D_2(2,1)=D_2(1,2);
D_2(3,1)=D_2(1,3);
D_2(2,3)=lambda(2);
D_2(3,2)=D_2(2,3);

% Store the tensor in a cell array
D{1}=D_1;
D{2}=D_2;

% Specify the material to be used a conductor (where gradients basis
% functions to be included)
% Also specify mechanical subdomain

matcond=[2];
subdom_mech=[1];
mur=mu(2)/mu(1);

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 0.01;  % Object size
shift=[0 0 0]; % Object shift

%-------------------------------------------------------------------------------
% Store material properties on probdata.matr substructure
%-------------------------------------------------------------------------------

probdata.matr.muz=muz;
probdata.matr.mur=mur;
probdata.matr.epz=epz;
probdata.matr.omega=omega;
probdata.matr.mu=mu;
probdata.matr.epl=epl;
probdata.matr.sigma=sigma;
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
probdata.matr.freq=freq;

%--------------------------------------------------------------------------
% Define prescribed constants
%--------------------------------------------------------------------------

B = 1;                                  % Magnitude of B_0
H0 = B/muz*[0;0;1];                     % Magnitude of H_0
rin(1) = 0.01;                          % Inner sphere radius
rin(2) = 1;                          % Outer sphere radius
p_0 = 1e9;                              % Pressure on the sphere


%-------------------------------------------------------------------------
% Blending Function Info (over write defaults)
%-------------------------------------------------------------------------
% Set geometry order
% 0 - linear
% 2 - quadratic
% 4 - Exact geometry for sphere

gorder =2;        % order (if gorder > 0, quadlin =2 required)

%Define parameters to check sphere approximation
sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once
svchk = 1;         % Check surface/volume as expecting a sphere


probdata.jb.gorder = gorder;
probdata.jb.rin = rin(1);
probdata.jb.sufv = sufv(1);
probdata.jb.gag = gag(1);
probdata.mesh.svchk=svchk;


% define boundary data

% current bctype
% bctype 5 Object (no bcs here) Interface boundary
% bctype 3 Neumann type (expected for this problem)
% bctype 2 Dirichlet (none here)



% Define constants needed for the solution of the harmonic EM problem




% Define arguments for exact solution and BC functions

arg.sigma=sigma(2);
arg.omega=omega;
arg.mu=mu(2)*muz;
arg.muz=muz;
arg.mur=mur;
arg.alpha=rin(1);
arg.B=B;
arg.H0=H0;
arg.muz=muz;
arg.rin=rin;
arg.nu=nu;
arg.E=E;
arg.p_a=p_0;
arg.p_b=p_0*1000;
arg.lambda=lambda(1);
arg.G=G(1);
arg.delta=delta;


% define the function handle for Dirichlet BC's
probdata.es.dirfun=@esproblemdir;
probdata.es.dirfunarg=arg;

% define the source terms
probdata.es.srcfunarg=[];
probdata.es.srcfun=@esproblemsrc;

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

%--------------------------------------------------------------------------

% Dirichlet Boundary Condition----------------------------------------------
function e=esproblemdir(x,y,z,index,arg,probstatic,probFlag)    

if probFlag==1                  % EM problem


if probstatic==1
   if index==2
e=[0;0;x*y];
   end

else



e=zeros(3,1);
if index==2
    % this would be non-zero for non zero type BC's
e=zeros(3,1);
else
    e(1)=complex(0,0);
    e(2)=complex(0,0);
    e(3)=complex(0,0);
end
end
elseif probFlag==2                      % Mechanical problem
    
  
    % Output vector
    e=[0;0;0];
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
   curle=[x;-y;0];

% Src term-----------------------------------------------------------------

function src=esproblemsrc(x,y,z,arg,domain,probstatic)               % Source term for es
src=[0;0;0];

function srcmech=esproblemsrcmech(x,y,z,arg,domain) % Source term for es
lambda=arg.lambda;
G=arg.G;
srcmech=[2*lambda+8*G;0;0];

% Analytical Solution Field------------------------------------------------
function [out]=esproblemexact(x,y,z,arg,probstatic,probFlag,omega)   

out=[0;0;x*y];

 
% Analytical Solution for B=curl A =mu H-----------------------------------
% Gives total field as output
function [curle,hf,H0]=esproblemexactcurl(x,y,z,arg,probstatic,omega)  

curle=[x;-y;0];
hf=zeros(3,1);
H0=zeros(3,1);

%--------------------------------------------------------------------------
% compute constants for the analytical solution
function [C,D]=myconstants(mu,omega,sigma,alpha,epsilon,mu0)

p=sigma*mu*omega;
v=sqrt(1i*p)*alpha;

% compute integrals

Imh=sqrt(2/pi/v)*cosh(v);


Iph=sqrt(2/pi/v)*sinh(v);

C=3*mu*v*alpha^(3/2)/( (mu-mu0)*v*Imh+(mu0*(1+v^2)-mu)*Iph);
D=alpha^3*((2*mu+mu0)*v*Imh-(mu0*(1+v^2)+2*mu)*Iph)/((mu-mu0)*v*Imh+(mu0*(1+v^2)-mu)*Iph);


function [efield,eddy]=electricfield(x,C,D,B,sigma,mu,omega,alpha)

[phi,thetamatlab,r] = cart2sph(x(1),x(2),x(3));
% thetamatlab is not the standard theta for spherical coordinates!
theta=pi/2-thetamatlab; % angle from z axis
% phi angle in x,y plane
% r distance from origin


p=sigma*mu*omega;
v=sqrt(1i*p)*r;
I32=(2/pi/v)^(1/2)*(cosh(v)-sinh(v)/v);

if r > alpha
    % no electric field outside the object
    Air=0;
    Aiphi=1/2*B*(r+D/r^2)*sin(theta);
    Aitheta=0;
else    
    Air=0;
    Aiphi=1/2*B*C/sqrt(r)*I32*sin(theta);
    Aitheta=0;
end

% Convert to Cartesian coordinates
Aix=sin(theta)*cos(phi)*Air+ cos(theta)*cos(phi)*Aitheta-sin(phi)*Aiphi;
Aiy=sin(theta)*sin(phi)*Air+ cos(theta)*sin(phi)*Aitheta+cos(phi)*Aiphi;
Aiz=cos(theta)*Air-sin(theta)*Aitheta;

% Convert to E-field
efield(1)=-1i*omega*Aix;
efield(2)=-1i*omega*Aiy;
efield(3)=-1i*omega*Aiz;


% Compute Eddy currents
if r <= alpha
    eddy=sigma*efield;
else
    eddy=0*efield;
end


% magnetic fields
function [bfield,hfield]=magneticfield(x,C,D,B,sigma,mu,omega,alpha,muz)

[phi,thetamatlab,r] = cart2sph(x(1),x(2),x(3));
% thetamatlab is not the standard theta for spherical coordinates!
theta=pi/2-thetamatlab; % angle from z axis
% phi angle in x,y plane
% r distance from origin


p=sigma*mu*omega;
k=sqrt(1i*p);
v=sqrt(1i*p)*r;
I32=(2/pi/v)^(1/2)*(cosh(v)-sinh(v)/v);
dI32r=sqrt(2/pi/k)*(-1/2)/r^(3/2)*(cosh(v)-sinh(v)/v)+(2/pi/v)^(1/2)*(sinh(v)*k-cosh(v)*k/v+sinh(v)/k/r^2);

if r >= alpha
    % B field outside the object
    Bitheta=-B*(1-D/2/r^3)*sin(theta);
    Biphi=0;
    Bir=B*(1+D/r^3)*cos(theta);
    Hitheta=Bitheta/muz;
    Hiphi=0;
    Hir=Bir/muz;

else    
    Bir=0.5*B*C/r^(3/2)*I32*2*cos(theta);
    Biphi=0;
    Bitheta=-1/r*(0.5*0.5*B*C/sqrt(r)*I32*sin(theta)+0.5*B*C*sqrt(r)*sin(theta)*dI32r);
    Hitheta=Bitheta/mu;
    Hiphi=0;
    Hir=Bir/mu;

end

% Convert to Cartesian coordinates
bfield(1)=sin(theta)*cos(phi)*Bir+ cos(theta)*cos(phi)*Bitheta-sin(phi)*Biphi;
bfield(2)=sin(theta)*sin(phi)*Bir+ cos(theta)*sin(phi)*Bitheta+cos(phi)*Biphi;
bfield(3)=cos(theta)*Bir-sin(theta)*Bitheta;


hfield(1)=sin(theta)*cos(phi)*Hir+ cos(theta)*cos(phi)*Hitheta-sin(phi)*Hiphi;
hfield(2)=sin(theta)*sin(phi)*Hir+ cos(theta)*sin(phi)*Hitheta+cos(phi)*Hiphi;
hfield(3)=cos(theta)*Hir-sin(theta)*Hitheta;

