function probdata=problem1ansys(pm)

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order

% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Polynomial Degree Info Order
order=pm; 


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
 vtkoption =1;
 
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


% %----------------------------------------------------
 %interpolation points for perturbed H
 
% %interpolation number
 N = 50;
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
 overp =[0.0,0,0.08];   % endpoint of the line
     
      point = [point; starp];   
     deltax = (overp(1)-starp(1))/(N2-1);
     deltay = (overp(2)-starp(2))/(N2-1);
     deltaz = (overp(3)-starp(3))/(N2-1);

     for i = 2:N2
         point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
     end 
     
     
 probdata.mesh.point = point;

%-------------------------------------------------------------------------

%job data------------------------------------------------------------------
job = 'gsphere'   % Job Name 
meshtype = 4;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
                    %ansys

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;
%------------------------------------------------------------------job data


probdata.sol.regopt =2 % 1- use regularisation with cg for gradent blocks
                       % 2- use regularisation with direct solve for gradient blocks
                       % 3- direct solve
%material data-------------------------------------------------------------
muz = 1.256637061435917e-06; % Mu_z
epz = 0;                     % Ep_z
omega =2*pi*133.5;           % Omega

% For each subdomain specify parameters
% Mu_r, Ep_r, Sigma, J
% Mat 1
mu(1) = 1;
epl(1) = 0;
sigma(1) = 0.1;
jsrc(1,1:3) = [0  0  0];

% Mat 2
mu(2) = 1.5
epl(2) = 0;
sigma(2) = 5.96e7;
jsrc(2,1:3) = [0  0  0];

% specify the material to be used a conductor conductors
% ie specify regions where gradients basis functions to be included
matcond=[2];

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 0.01;  % Object size
shift=[0 0 0]; % Object shift

probdata.matr.muz=muz;
probdata.matr.epz=epz;
probdata.matr.omega=omega;
probdata.matr.mu=mu;
probdata.matr.epl=epl;
probdata.matr.sigma=sigma;
probdata.matr.jsrc=jsrc;
probdata.matr.delta=delta;
probdata.matr.shift=shift;
probdata.matr.matcond=matcond;

%-------------------------------------------------------------------------
% Blending Function Info (over write defaults)
% gorder =4;        % order (if gorder > 0, quadlin =2 required)
% g1 = 2;           % use exact geometry for a sphere
%  rin(1) = 1;        % raduis
%  sufv(1)= 5;        % surfaceval
%  gag = 2;           % Goagain  % 1- go again 2- once
%  svchk = 1;         % Check surface/volume as expecting a sphere
% % 
% % 
% % probdata.jb.gorder = gorder;
% % probdata.jb.g1 = g1;
%  probdata.jb.rin = rin(1);
%  probdata.jb.sufv = sufv(1);
%  probdata.jb.gag = gag(1);
%  probdata.mesh.svchk=svchk;

% define boundary data

% current bctype
% bctype 5 Object (no bcs here)
% bctype 3 Neumann type (expected for this problem)
% bctype 2 Dirichlet (none here)

B=muz; % magnitude of B_z
[C,D]=myconstants(mu(2)*muz,omega,sigma(2),delta,0,muz);
arg.sigma=sigma(2);
arg.omega=omega;
arg.mu=mu(2)*muz;
arg.muz=muz;
arg.alpha=delta;
arg.C=C;
arg.D=D;
arg.B=B;

% define the function handle for Dirichlet BC's
probdata.es.dirfun=@esproblemdir;
probdata.es.dirfunarg=arg;

% define the source terms
probdata.es.srcfunarg=[];
probdata.es.srcfun=@esproblemsrc;

% define the function handle for Neumann BC's
H0=B/muz*[0;0;1];
arg.H0=H0;
arg.muz=muz;
%arg.mtensor=mtensor
probdata.es.neufun=@esproblemneu;
probdata.es.neufunarg=arg;

% exact solution
probdata.es.exactfun=@esproblemexact; 
probdata.es.exactfunarg=arg;

probdata.es.exactcurlfun=@esproblemexactcurl;
probdata.es.exactcurlfunarg=arg;

%--------------------------------------------------------------------------

% Dirichlet Boundary Conditio----------------------------------------------
function e=esproblemdir(x,y,z,index,arg)       % DBC for es (4*pi*10e-7)^(0.5)*
sigma=arg.sigma;
omega=arg.omega;
mu=arg.mu;
muz=arg.muz;
alpha=arg.alpha;
C=arg.C;
D=arg.D;
B=arg.B;


[efield,eddy]=electricfield([x,y,z],C,D,B,sigma,mu,omega,alpha);
e=zeros(3,1);
if index==2
    % this would be non-zero for non zero type BC's
     e(1)=efield(1)./(-1i*omega);
     e(2)=efield(2)./(-1i*omega);
     e(3)=efield(3)./(-1i*omega);
else
    e(1)=complex(0,0);
    e(2)=complex(0,0);
    e(3)=complex(0,0);
end

% Neumann Boundary Condition-----------------------------------------------

function curle=esproblemneu(x,y,z,index,arg)
if index==3
sigma=arg.sigma;
omega=arg.omega;
mu=arg.mu;
muz=arg.muz;
alpha=arg.alpha;
C=arg.C;
D=arg.D;
B=arg.B;

[bfield,hfield]=magneticfield([x,y,z],C,D,B,sigma,mu,omega,alpha,muz);
 curle=zeros(3,1);
 curle(1)=bfield(1);
 curle(2)=bfield(2);
 curle(3)=bfield(3);
else
    curle=zeros(3,1);
end    

% Src term-----------------------------------------------------------------

function src=esproblemsrc(x,y,z,dum,arg)               % Source term for es
src=[0;0;0];

% Analytical Solution Field------------------------------------------------
function [out]=esproblemexact(x,y,z,arg)           
sigma=arg.sigma;
omega=arg.omega;
mu=arg.mu;
muz=arg.muz;
alpha=arg.alpha;
C=arg.C;
D=arg.D;
B=arg.B;


[efield,eddy]=electricfield([x,y,z],C,D,B,sigma,mu,omega,alpha);
out=(efield.')./(-1i*omega);

 
% Analytical Solution for B=curl A =mu H-----------------------------------
% Gives total field as output
function [curle,hf,H0]=esproblemexactcurl(x,y,z,arg)         

sigma=arg.sigma;
omega=arg.omega;
mu=arg.mu;
muz=arg.muz;
alpha=arg.alpha;
C=arg.C;
D=arg.D;
B=arg.B;

[bfield,hfield]=magneticfield([x,y,z],C,D,B,sigma,mu,omega,alpha,muz);
curle=bfield.';
hf=hfield.';
H0=B/muz*[0 0 1].';

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


