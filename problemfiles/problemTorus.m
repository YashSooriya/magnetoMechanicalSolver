function probdata=problemTorus(pm)

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
 erroroption =0; % Exact E field not entered!
 
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


% %----------------------------------------------------
 %interpolation points for perturbed H
 
% %interpolation number
 N = 100;
 probdata.mesh.N = N;
 N1=20;
 N2=N-N1;

 starp = [0.0,0,0.00];   % startpoin of the line
 overp =[0.009,0.009,0.009];   % endpoint of the line
      point = starp;   
     deltax = (overp(1)-starp(1))/(N1-1);
     deltay = (overp(2)-starp(2))/(N1-1);
     deltaz = (overp(3)-starp(3))/(N1-1);

     for i = 2:N1
         point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
     end 
 starp = [0,0,0.011];   % startpoin of the line
 overp =[0.25,0.25,0.25];   % endpoint of the line
     
      point = [point; starp];   
     deltax = (overp(1)-starp(1))/(N2-1);
     deltay = (overp(2)-starp(2))/(N2-1);
     deltaz = (overp(3)-starp(3))/(N2-1);

     for i = 2:N2
         point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
     end 
     
     
 probdata.mesh.point = point;


%job data------------------------------------------------------------------
%job data------------------------------------------------------------------
job='torus2';

% job='sphere1';

meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;

%------------------------------------------------------------------job data


probdata.sol.regopt =2; % 1- use regularisation with cg for gradent blocks
                       % 2- use regularisation with direct solve for gradient blocks
                       % 3- direct solve
%material data-------------------------------------------------------------
nmat = 2;         % Number of materieals
muz = 1.256637061435917e-06; % Mu_z
epz = 0;                     % Ep_z
omega = 133.5;               % Omega

% Mu_r, Ep_r, Sigma, J
% Mat 1
mu(1) = 1;
epl(1) = 0;
sigma(1) =0.1;
jsrc(1,1:3) = [0  0  0];
% Mat 2
mu(2) = 1;
epl(2) = 0;
sigma(2) = 5.96e7;
jsrc(2,1:3) = [0  0  0];

% specify the material to be used a conductor conductors
% ie specify regions where gradients basis functions to be included
matcond=[2];

% The object has been generated 100 times bigger in order to avoid a
% problem with NetGen. Now we re-escale the mesh with the factor delta.
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
gorder =4;        % order (if gorder > 0, quadlin =2 required)
g1 = 4;           

sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once
svchk = 0;         % Check surface/volume as expecting a sphere
rin(1)=0.01;
rin(2)=0.02;
probdata.jb.rin = rin(1);


probdata.jb.gorder = gorder;
probdata.jb.g1 = g1;

probdata.jb.sufv = sufv(1);
probdata.jb.gag = gag(1);
probdata.mesh.svchk=svchk;


% define boundary data

% current bctype
% bctype 5 Object (no bcs here)
% bctype 3 Neumann type (expected for this problem)
% bctype 2 Dirichlet (none here)
arg=[];

% define the function handle for Dirichlet BC's
probdata.es.dirfun=@esproblemdir;
probdata.es.dirfunarg=arg;

% define the source terms
probdata.es.srcfunarg=[];
probdata.es.srcfun=@esproblemsrc;

% define the function handle for Neumann BC's

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
e=zeros(3,1);


% Neumann Boundary Condition-----------------------------------------------

function curle=esproblemneu(x,y,z,index,arg)
curle=zeros(3,1);
curle(1)=1;

% Src term-----------------------------------------------------------------

function src=esproblemsrc(x,y,z,dum,arg,domain)               % Source term for es
src=[0;0;0];

