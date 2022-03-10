function probdata=problemTEAM7(pm)

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
 plotoption = 0;
 
% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 0;
 
% Calculate the H(curl) norm of the error
% % erroroption = 0 - do not calculate error
% % erroroption = 1 - calculate error
 erroroption = 0; % Exact E field not entered!
 
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



%job data------------------------------------------------------------------
%job data------------------------------------------------------------------

job='team7_2';

meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.order=order;
probdata.jb.name='Team7_coil';


 
% %interpolation number
 Nb = 100;
 Ne = 100;
 probdata.mesh.Nb = Nb;
 probdata.mesh.Ne = Ne;


 starp = [0.0001,0.072,0.034];   % startpoin of the line
 overp =[0.294,0.072,0.034];   % endpoint of the line
      point = starp;   
     deltax = (overp(1)-starp(1))/(Nb-1);
     deltay = (overp(2)-starp(2))/(Nb-1);
     deltaz = (overp(3)-starp(3))/(Nb-1);

     for i = 2:Nb
         point = [point;starp(1)+deltax*(i-1) starp(2)+deltay*(i-1) starp(3)+deltaz*(i-1)];    %coordinate of interpolation points
     end 
  % ============================================================= point 2 =============================================================   
     N1=10;
     N2=90;
     
      starp2 = [0.0001,0.072,0.0189];   % startpoin of the line
 overp2 =[0.0179,0.072,0.0189];   % endpoint of the line
      point2 = starp2;   
     deltax2 = (overp2(1)-starp2(1))/(N1-1);
     deltay2 = (overp2(2)-starp2(2))/(N1-1);
     deltaz2 = (overp2(3)-starp2(3))/(N1-1);

     for i = 2:N1
         point2 = [point2;starp2(1)+deltax2*(i-1) starp2(2)+deltay2*(i-1) starp2(3)+deltaz2*(i-1)];    %coordinate of interpolation points
     end 
           starp2 = [0.1261,0.072,0.0189];   % startpoin of the line
 overp2 =[0.293,0.072,0.0189];   % endpoint of the line
      point2 = [point2;starp2];   
     deltax2 = (overp2(1)-starp2(1))/(N2-1);
     deltay2 = (overp2(2)-starp2(2))/(N2-1);
     deltaz2 = (overp2(3)-starp2(3))/(N2-1);

     for i = 2:N2
         point2 = [point2;starp2(1)+deltax2*(i-1) starp2(2)+deltay2*(i-1) starp2(3)+deltaz2*(i-1)];    %coordinate of interpolation points
     end 
%=========================================== point3 ======================================================================================
         starp3 = [0.0001,0.072,0.0001];   % startpoin of the line
 overp3 =[0.0179,0.072,0.0001];   % endpoint of the line
      point3 = starp3;   
     deltax3 = (overp3(1)-starp3(1))/(N1-1);
     deltay3 = (overp3(2)-starp3(2))/(N1-1);
     deltaz3 = (overp3(3)-starp3(3))/(N1-1);

     for i = 2:N1
         point3 = [point3;starp3(1)+deltax3*(i-1) starp3(2)+deltay3*(i-1) starp3(3)+deltaz3*(i-1)];    %coordinate of interpolation points
     end 
           starp3 = [0.1261,0.072,0.0001];   % startpoin of the line
 overp3 =[0.293,0.072,0.0001];   % endpoint of the line
      point3 = [point3;starp3];   
     deltax3 = (overp3(1)-starp3(1))/(N2-1);
     deltay3 = (overp3(2)-starp3(2))/(N2-1);
     deltaz3 = (overp3(3)-starp3(3))/(N2-1);

     for i = 2:N2
         point3 = [point3;starp3(1)+deltax3*(i-1) starp3(2)+deltay3*(i-1) starp3(3)+deltaz3*(i-1)];    %coordinate of interpolation points
     end 
     
 % ============================================================ point 4 =============================================================
      starp4 = [0.0001,0.144,0.034];   % startpoin of the line
 overp4 =[0.294,0.144,0.034];   % endpoint of the line
      point4 = starp4;   
     deltax4 = (overp4(1)-starp4(1))/(Nb-1);
     deltay4 = (overp4(2)-starp4(2))/(Nb-1);
     deltaz4 = (overp4(3)-starp4(3))/(Nb-1);

     for i = 2:Nb
         point4 = [point4;starp4(1)+deltax4*(i-1) starp4(2)+deltay4*(i-1) starp4(3)+deltaz4*(i-1)];    %coordinate of interpolation points
     end 

     
     
 probdata.mesh.point = point;
 probdata.mesh.point2 = point2;
 probdata.mesh.point3 = point3;
 probdata.mesh.point4 = point4;

%--------------------------------------------------------------------------


probdata.sol.regopt =1; % 1- use regularisation with cg for gradent blocks
                       % 2- use regularisation with direct solve for gradient blocks
                       % 3- direct solve
%material data-------------------------------------------------------------
%nmat = 3;         % Number of materieals
muz = 1.256637061435917e-06; % Mu_0
epz = 0;                     % Ep_0
omega = 100*3.141592;               % Omega

% Mu_r, Ep_r, Sigma, J
% Mat conductor
mu(1:19) = 1;
epl(1:19) = 0;
sigma(1) = 3.526e7;
jsrc(1,1:3) = [0  0  0];
% Mat coil
for imat=2:9
sigma(imat) = 0.1;
jsrc(imat,1:3) = [0  0  0];
end
% Air
for imat=10:19
    sigma(imat)=0.1;
    jsrc(imat,1:3) = [0 0 0];
end


% specify the material to be used a conductor conductors
% ie specify regions where gradients basis functions to be included
matcond=[1 2:9];



probdata.matr.muz=muz;
probdata.matr.epz=epz;
probdata.matr.omega=omega;
probdata.matr.mu=mu;
probdata.matr.epl=epl;
probdata.matr.sigma=sigma;
probdata.matr.jsrc=jsrc;
probdata.matr.matcond=matcond;


%-------------------------------------------------------------------------
% Blending Function Info (over write defaults)
gorder =4;        % order (if gorder > 0, quadlin =2 required)
g1 = 4;           

sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once


probdata.jb.gorder = gorder;
probdata.jb.g1 = g1;
probdata.jb.sufv = sufv(1);
probdata.jb.gag = gag(1);

delta=1;
shift=[0 0 0]; % Object shift
probdata.matr.delta=delta;
probdata.matr.shift=shift;

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


%arg.mtensor=mtensor
probdata.es.neufun=@esproblemneu;
probdata.es.neufunarg=arg;



%--------------------------------------------------------------------------

% Dirichlet Boundary Conditio----------------------------------------------
function e=esproblemdir(x,y,z,index,arg)       % DBC for es (4*pi*10e-7)^(0.5)*
e=zeros(3,1);


% Neumann Boundary Condition-----------------------------------------------

function curle=esproblemneu(x,y,z,index,arg)
curle=zeros(3,1);


% Src term-----------------------------------------------------------------
function src=esproblemsrc(x,y,z,dum,arg,domain)              % Source term for es
Jmag=1.0968e6*1.256637061435917e-06;
N=length(x);
for i=1:N
if domain==3 
    x(i)=x(i)-0.194-(0.075-0.025/sqrt(2));
    y(i)=y(i)-0.1-(0.075-0.025/sqrt(2));
    theta=atan2(y(i),x(i));
src=Jmag*[-sin(theta);cos(theta);0];
elseif domain==9
    x(i)=x(i)-0.194-(0.075-0.025/sqrt(2));
    y(i)=y(i)-0.1+(0.075-0.025/sqrt(2));
    theta=atan2(y(i),x(i));
    src=Jmag*[-sin(theta);cos(theta);0];
elseif domain==5 
    y(i)=y(i)-0.1-(0.075-0.025/sqrt(2));
    x(i)=x(i)-0.194+(0.075-0.025/sqrt(2));
    theta=atan2(y(i),x(i));
    src=Jmag*[-sin(theta);cos(theta);0];
elseif domain==7
    x(i)=x(i)-0.194+(0.075-0.025/sqrt(2));
    y(i)=y(i)-0.1+(0.075-0.025/sqrt(2));
    theta=atan2(y(i),x(i));
    src=Jmag*[-sin(theta);cos(theta);0];
elseif domain==2 
    src=+Jmag*[0;1;0];
elseif domain==6
    src=-Jmag*[0;1;0];
elseif domain==4
    src=-Jmag*[1;0;0];
elseif domain==8
    src=Jmag*[1;0;0];
else
    src=[0;0;0];
end
end


%--------------------------------------------------------------------------








