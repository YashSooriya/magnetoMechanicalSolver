function probdata=problemCylinder(pm)

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order

% include standard options (include BC, src definition etc)
probdata=[];
probdata=standoptions(probdata);

% Choose the BC type for the outer boundary
probdata.bctypeouter=2;   % 2 Dirichlet outer BC
                          % 3 Neumann outer BC

% Polynomial Degree Info Order
order=pm; 

%--------------------------------------------------------------------------
% Define postprocessing switchs
%--------------------------------------------------------------------------

% % plot the perturbed solution on a line
% % plotoption = 0 - do not plot
% % plotoption = 1 - plot
 plotoption = 0;
 
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



%--------------------------------------------------------------------------
% Job data
%--------------------------------------------------------------------------

job=sprintf('cylinder_lin_1');   % Job Name 
%job=sprintf('cylinder_lin_pref');


% job='sphere1';3
meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
                   % 4= ansys, 5 = Opera
probdata.jb.name=job;
%probdata.jb.name=sprintf('cylinder_lin_3');
probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;

%--------------------------------------------------------------------------
% Define solver option
%--------------------------------------------------------------------------


probdata.sol.regopt =2;     % 1- use regularisation with cg for gradent blocks
                            % 2- use regularisation with direct solve for gradient blocks
                            % 3- direct solve
                            
TOL_GMRES=1e-5;
probdata.TOL_GMRES=TOL_GMRES;
                            
%--------------------------------------------------------------------------
% Material properties
%--------------------------------------------------------------------------

nmat = 1;                            % Number of materials
muz = 1.256637061435917e-06;         % Mu_z (Free sopace permeability)
epz = 0;                             % Ep_z
omega = 50;

% Mu_r, Ep_r, Sigma, J

% Mat 1
mu(1) = 1;                           % Permeability (relative)
epl(1) = 0;                          % Permittivity (relative)
sigma(1) = 1e-3;                     % Conductivity
regTerm=1e-3;                        % Regularization term
jsrc(1,1:3) = [0  0  0];
% Mechanical properties
rho(1) = 7800;                       % Density
E(1) = 210e9;                        % Young's modulus
%E(1)=1000;
nu(1) = 0.33;                        % Poisson's ratio



% Define Lame constants
lambda   = (E.*nu)./((1+nu).*(1-2.*nu));
G        = E./(2.*(1+nu)); 

% Define the elasticity tensor(s). We'll have as many as different material

% For this problem we have just one
D_1=zeros(6);
D_1(1,1)=lambda+2*G;
D_1(2,2)=lambda+2*G;
D_1(3,3)=lambda+2*G;
D_1(4,4)=G;
D_1(5,5)=G;
D_1(6,6)=G;
D_1(1,2)=lambda;
D_1(1,3)=lambda;
D_1(2,1)=D_1(1,2);
D_1(3,1)=D_1(1,3);
D_1(2,3)=lambda;
D_1(3,2)=D_1(2,3);

% Store the tensor in a cell array
D{1}=D_1;


% Specify the material to be used a conductor (where gradients basis
% functions to be included)
% Also specify mechanical subdomain

matcond=[1];
subdom_mech=[1];
mur=mu(1);

% In this case the mesh is for a unit sized object, it must be scaled (and
% repositioned)
delta = 1;  % Object size
shift=[0 0 2.5]; % Object shift
matMech=1;
NmechBodies=1;
matAir=2;

%-------------------------------------------------------------------------------
% Store material properties on probdata.matr substructure
%-------------------------------------------------------------------------------

probdata.matr.muz=muz;
probdata.matr.mur=mur;
probdata.matr.epz=epz;
probdata.matr.mu=mu;
probdata.matr.epl=epl;
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
probdata.matMech=matMech;
probdata.NmechBodies=NmechBodies;
probdata.matAir=matAir;

%--------------------------------------------------------------------------
% Define prescribed constants
%--------------------------------------------------------------------------

B = 1;                                  % Magnitude of B_0
H0 = B/muz*[0;0;1];                     % Magnitude of H_0
rin(1) = 1;                          % Inner sphere radius
rin(2) = 2;                          % Outer sphere radius
p_0 = 1e4;                              % Pressure on the sphere


%-------------------------------------------------------------------------
% Blending Function Info (over write defaults)
%-------------------------------------------------------------------------
gorder =0;        % order (if gorder > 0, quadlin =2 required)
g1 = 1;           % use exact geometry for a sphere

sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once
svchk = 1;         % Check surface/volume as expecting a sphere


probdata.jb.gorder = gorder;
probdata.jb.g1 = g1;
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

[C,D]=myconstants(mu(1)*muz,omega,sigma(1),delta,0,muz);


% Define arguments for exact solution and BC functions

arg.sigma=sigma(1);
arg.mu=mu(1)*muz;
arg.muz=muz;
arg.mur=mur;
arg.alpha=delta;
arg.C=C;
arg.D=D;
arg.B=B;
arg.H0=H0;
arg.muz=muz;
arg.rin=delta*rin;
arg.nu=nu;
arg.E=E;
arg.p_a=p_0;
arg.p_b=p_0*1000;


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
function e=esproblemdir(x,y,z,index,arg,probstatic,probFlag,omega)    

if probFlag==1                  % EM problem
sigma=arg.sigma;

omega=arg.omega;
mu=arg.mu;
muz=arg.muz;
mur=arg.mur;
alpha=arg.alpha;
C=arg.C;
D=arg.D;

B=arg.B;

if probstatic==1
    if index==2
    A=zeros(3,1);
    r=sqrt(x^2+y^2+z^2);
    theta=acos(z/r);
    phi=atan2(y,x);
    if r>alpha
    Aphi=0.5*B*(r-2*(muz-mur*muz)/(2*muz+mur*muz)*alpha^3/r^2)*sin(theta);
  
    else
    Aphi=3*B*r*sin(theta)/(2*(1+2/mur));
    end
    A(1)=-sin(phi)*Aphi;
    A(2)=cos(phi)*Aphi;
    A(3)=0;
    e=A;
    else
         e(1)=complex(0,0);
    e(2)=complex(0,0);
    e(3)=complex(0,0);
    end
else


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
end
elseif probFlag==2                      % Mechanical problem
    
    % Extract problem constants
    E=arg.E(1);
    nu=arg.nu(1);
    b=arg.rin(2);
    a=arg.rin(1);
    p_a=arg.p_a;
    p_b=arg.p_b;
    
    % Define spherical coordinates
    R=sqrt(x^2+y^2);

    phi=atan2(y,x);
    
    % Analytical solution in cylindrical coordinates
    u_r=(1+nu)*a^2*b^2/(E*(b^2-a^2))*((p_a-p_b)/R+(1-2*nu)*(p_a*a^2-p_b*b^2)*R/(a^2*b^2))-nu*R*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+2*nu^2*R*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2));
    u_z=z*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))-2*nu*z*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2));
    % Analytical solution in Cartesian coordinates
    u_x=cos(phi)*u_r;
    u_y=sin(phi)*u_r;
  
    % Output vector
    e=[u_x;u_y;u_z];
end

function grade=gradExact(x,y,z,arg)

  E=arg.E(1);
    nu=arg.nu(1);
    b=arg.rin(2);
    a=arg.rin(1);
    p_a=arg.p_a;
    p_b=arg.p_b;
    R=sqrt(x^2+y^2);
    phi=atan2(y,x);
  
   
%     a_11=(-sin(atan2(y,x))*(-(imag(x) + real(y))/((imag(x) + real(y))^2 + (imag(y) - real(x))^2)))*(1+nu)*a^2*b^2/(E*(b^2-a^2))*((p_a-p_b)/sqrt(x^2+y^2)+(1-2*nu)*(p_a*a^2-p_b*b^2)*sqrt(x^2+y^2)/(a^2*b^2))-nu*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+2*nu^2*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+(-(imag(y) - real(x))/abs(x + y*1i))*(nu*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((x*(p_a - p_b))/(x^2 + y^2)^(3/2) + (x*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2));
%     a_12=(-sin(atan2(y,x))*(-(imag(y) - real(x))/((imag(x) + real(y))^2 + (imag(y) - real(x))^2)))*(1+nu)*a^2*b^2/(E*(b^2-a^2))*((p_a-p_b)/sqrt(x^2+y^2)+(1-2*nu)*(p_a*a^2-p_b*b^2)*sqrt(x^2+y^2)/(a^2*b^2))-nu*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+2*nu^2*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+(imag(x) + real(y))/abs(x + y*1i)*(((nu*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((y*(p_a - p_b))/(x^2 + y^2)^(3/2) + (y*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2))));
%     a_21=(cos(atan2(y,x))*(-(imag(x) + real(y))/((imag(x) + real(y))^2 + (imag(y) - real(x))^2)))*(1+nu)*a^2*b^2/(E*(b^2-a^2))*((p_a-p_b)/sqrt(x^2+y^2)+(1-2*nu)*(p_a*a^2-p_b*b^2)*sqrt(x^2+y^2)/(a^2*b^2))-nu*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+2*nu^2*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+((imag(x) + real(y))/abs(x + y*1i))*((nu*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((x*(p_a - p_b))/(x^2 + y^2)^(3/2) + (x*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)));
%     a_22=(cos(atan2(y,x))*(-(imag(y) - real(x))/((imag(x) + real(y))^2 + (imag(y) - real(x))^2)))*((1+nu)*a^2*b^2/(E*(b^2-a^2))*((p_a-p_b)/sqrt(x^2+y^2)+(1-2*nu)*(p_a*a^2-p_b*b^2)*sqrt(x^2+y^2)/(a^2*b^2))-nu*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+2*nu^2*sqrt(x^2+y^2)*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2)))+((imag(x) + real(y))/abs(x + y*1i))*((nu*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((y*(p_a - p_b))/(x^2 + y^2)^(3/2) + (y*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)));
%     a_11=((nu*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((x*(p_a - p_b))/(x^2 + y^2)^(3/2) + (x*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)))/(y^2/x^2 + 1)^(1/2) - (y^2*((2*nu^2*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (nu*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*((p_a - p_b)/(x^2 + y^2)^(1/2) - ((x^2 + y^2)^(1/2)*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2))*(nu + 1))/(E*(a^2 - b^2))))/(x^3*(y^2/x^2 + 1)^(3/2));
%   a_12=((nu*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((y*(p_a - p_b))/(x^2 + y^2)^(3/2) + (y*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)))/(y^2/x^2 + 1)^(1/2) + (y*((2*nu^2*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (nu*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*((p_a - p_b)/(x^2 + y^2)^(1/2) - ((x^2 + y^2)^(1/2)*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2))*(nu + 1))/(E*(a^2 - b^2))))/(x^2*(y^2/x^2 + 1)^(3/2));
%   a_13=0;
%   a_21=(y*((2*nu^2*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (nu*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*((p_a - p_b)/(x^2 + y^2)^(1/2) - ((x^2 + y^2)^(1/2)*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2))*(nu + 1))/(E*(a^2 - b^2))))/(x^2*(y^2/x^2 + 1)^(1/2)) - (y^3*((2*nu^2*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (nu*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*((p_a - p_b)/(x^2 + y^2)^(1/2) - ((x^2 + y^2)^(1/2)*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2))*(nu + 1))/(E*(a^2 - b^2))))/(x^4*(y^2/x^2 + 1)^(3/2)) + (y*((nu*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((x*(p_a - p_b))/(x^2 + y^2)^(3/2) + (x*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*x*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2))))/(x*(y^2/x^2 + 1)^(1/2));
%   a_22=(y^2*((2*nu^2*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (nu*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*((p_a - p_b)/(x^2 + y^2)^(1/2) - ((x^2 + y^2)^(1/2)*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2))*(nu + 1))/(E*(a^2 - b^2))))/(x^3*(y^2/x^2 + 1)^(3/2)) - ((2*nu^2*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (nu*(x^2 + y^2)^(1/2)*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*((p_a - p_b)/(x^2 + y^2)^(1/2) - ((x^2 + y^2)^(1/2)*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2))*(nu + 1))/(E*(a^2 - b^2)))/(x*(y^2/x^2 + 1)^(1/2)) + (y*((nu*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2)) + (a^2*b^2*((y*(p_a - p_b))/(x^2 + y^2)^(3/2) + (y*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2*(x^2 + y^2)^(1/2)))*(nu + 1))/(E*(a^2 - b^2)) - (2*nu^2*y*(p_a*a^2 - p_b*b^2))/(E*(x^2 + y^2)^(1/2)*(a^2 - b^2))))/(x*(y^2/x^2 + 1)^(1/2));
%   a_23=0;
%   a_31=0;
%   a_32=0;
%   a_33=(2*nu*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (p_a*a^2 - p_b*b^2)/(E*(a^2 - b^2));

% define the gradient in cylindrical coord
a_11=(nu*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (2*nu^2*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*((p_a - p_b)/R^2 + ((2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2))*(nu + 1))/(E*(a^2 - b^2));
a_12=0;
a_13=0;
a_21=0;
a_22=-((2*R*nu^2*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (R*nu*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) + (a^2*b^2*(nu + 1)*((p_a - p_b)/R - (R*(2*nu - 1)*(p_a*a^2 - p_b*b^2))/(a^2*b^2)))/(E*(a^2 - b^2)))/R;
a_23=0;
a_31=0;
a_32=0;
a_33=(2*nu*(p_a*a^2 - p_b*b^2))/(E*(a^2 - b^2)) - (p_a*a^2 - p_b*b^2)/(E*(a^2 - b^2));

  gradecyl=[a_11 a_12 a_13;a_21 a_22 a_23;a_31 a_32 a_33];
% Transform this tensor to cartesian coord
grade=[cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0;0 0 1]*gradecyl*[cos(phi) sin(phi) 0;-sin(phi) cos(phi) 0;0 0 1];


 


% Neumann Boundary Condition-----------------------------------------------
probstatic=1;
function curle=esproblemneu(x,y,z,index,arg,probstatic,omega)
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

function src=esproblemsrc(x,y,z,dum,arg,domain)               % Source term for es
src=[0;0;0];

function srcmech=esproblemsrcmech(x,y,z,dum,arg,domain) % Source term for es
srcmech=[0;0;0];

% Analytical Solution Field------------------------------------------------
function [out]=esproblemexact(x,y,z,arg,probstatic,probFlag,omega)   

if probFlag==1              % EM problem
sigma=arg.sigma;


omega=arg.omega;
mu=arg.mu;
mur=arg.mur;
muz=arg.muz;
alpha=arg.alpha;
C=arg.C;
D=arg.D;
B=arg.B;


if probstatic==1
    A=zeros(3,1);
    r=sqrt(x^2+y^2+z^2);
    theta=acos(z/r);
    phi=atan2(y,x);
    if r>alpha
          Aphi=0.5*B*(r-2*(muz-mur*muz)/(2*muz+mur*muz)*alpha^3/r^2)*sin(theta);
    else
          Aphi=3*B*r*sin(theta)/(2*(1+2/mur));
    end
    A(1)=-sin(phi)*Aphi;
    A(2)=cos(phi)*Aphi;
    A(3)=0;
    out=A.';
else
[efield,eddy]=electricfield([x,y,z],C,D,B,sigma,mu,omega,alpha);
out=(efield.')./(-1i*omega);
end
elseif probFlag==2                  % Mechanical problem
    
    % Extract problem constants
    E=arg.E(1);
    b=arg.rin(2);
    a=arg.rin(1);
    p_a=arg.p_a;
    p_b=arg.p_b;
    nu=arg.nu(1);
    
     % Define spherical coordinates
    R=sqrt(x^2+y^2);

    phi=atan2(y,x);
    
    % Analytical solution in spherical coordinates
    u_r=(1+nu)*a^2*b^2/(E*(b^2-a^2))*((p_a-p_b)/R+(1-2*nu)*(p_a*a^2-p_b*b^2)*R/(a^2*b^2))-nu*R*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))+2*nu^2*R*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2));
    u_z=z*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2))-2*nu*z*(p_a*a^2-p_b*b^2)/(E*(b^2-a^2));
    % Analytical solution in Cartesian coordinates
    u_x=cos(phi)*u_r;
    u_y=sin(phi)*u_r;
  
    % Output vector
    out=[u_x;u_y;u_z];
end

 
% Analytical Solution for B=curl A =mu H-----------------------------------
% Gives total field as output
function [curle,hf,H0]=esproblemexactcurl(x,y,z,arg,probstatic,omega)  

sigma=arg.sigma;

omega=arg.omega;
mu=arg.mu;
mur=arg.mur;
muz=arg.muz;
alpha=arg.alpha;
C=arg.C;
D=arg.D;
B=arg.B;

if probstatic==1
    r=sqrt(x^2+y^2+z^2);
    phi=atan2(y,x);
    theta=acos(z/r);

    if r>alpha
    Br=B*(1-2*(muz-mur*muz)/(2*muz+mur*muz)*alpha^3/r^3)*cos(theta);
    Btheta=-B*sin(theta)*(1-(muz-mur*muz)/(2*muz+mur*muz)*alpha^3*log(r)/r);
    Hr=Br/muz;
    Htheta=Btheta/muz;
    
    else
    Br=3*B*cos(theta)/(1+2/mur);
    Btheta=-3*B*sin(theta)/(1+2/mur);
    Hr=Br/(mur*muz);
    Htheta=Btheta/(mur*muz);
    end
    Bfield(1)=sin(theta)*cos(phi)*Br+cos(theta)*cos(phi)*Btheta;
    Bfield(2)=sin(theta)*sin(phi)*Br+cos(theta)*sin(phi)*Btheta;
    Bfield(3)=cos(theta)*Br-sin(theta)*Btheta;
    H(1)=sin(theta)*cos(phi)*Hr+cos(theta)*cos(phi)*Htheta;
    H(2)=sin(theta)*sin(phi)*Hr+cos(theta)*sin(phi)*Htheta;
    H(3)=cos(theta)*Hr-sin(theta)*Htheta;
    curle=Bfield.';
    hf=H.';
H0=B/muz*[0 0 1].';
else
[bfield,hfield]=magneticfield([x,y,z],C,D,B,sigma,mu,omega,alpha,muz);
curle=bfield.';
hf=hfield.';
H0=B/muz*[0 0 1].';
end

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

