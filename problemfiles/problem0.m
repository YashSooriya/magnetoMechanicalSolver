function probdata=problem0(pm)

% Set the polynomial degree of the elements
% order = 0, 1, 2, 3
% nb order=0 refes to the standard linear/bilinear hat functions
% start adapting from this order

% Polynomial Degree Info Order
order=pm; 

% % plot the solution on a line
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
 erroroption = 1;
 
% % output the VTK file
% % vtkoption = 0 - do not output
% % vtkoption = 1 - output
 vtkoption =1;

% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 1; 

% % plot the  H_alpha-H_0 - D2G PT H_0  on a line
% % darkooption = 0 - do not plot
% % darkooption = 1 - plot
 darkooption = 1; 
 
probdata.mesh.plotoption=plotoption;
probdata.mesh.erroroption=erroroption;
probdata.mesh.vtkoption=vtkoption;
probdata.mesh.voltoption=voltoption;
probdata.mesh.darkooption=darkooption;



%job data------------------------------------------------------------------
job = 'sphere11_dir' %'coilwithsphere6neu' %'sphere6_2' %
meshtype = 3;      % Mesh type 1= FLITE , 2=NG (old style), 3=NG (new style)
matflg = 1;        % Mat Flag from file
ppus = 1 ;         % Parrellel Processors (not used)
% Polynomial Degree Info
ptype = 2;         % (Type)
% Blending Function Info
gorder = 0;        % order (if gorder > 0, quadlin =2 required)
g1 = 4;            % type (g1=4 recomended)
rin(1) = 1;        % raduis (used only if g1=2)
sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once
quadlin=1;         % quadlin % 1 linear , 2% quadratic geometry in netgen 
                   % mesh data file

probdata.jb.job = job;
probdata.jb.meshtype = meshtype;
probdata.jb.matflg=matflg;
probdata.jb.ppus = ppus;
probdata.jb.ptype = ptype; 
probdata.jb.order=order;
probdata.jb.gorder = gorder;
probdata.jb.g1 = g1;
probdata.jb.rin = rin(1);
probdata.jb.sufv = sufv(1);
probdata.jb.gag = gag(1);
probdata.jb.nrhs=1;

probdata.jb.quadlin = quadlin;
%------------------------------------------------------------------job data


%material data-------------------------------------------------------------
nmat = 2;                    % Number of materieals
muz = 1;                     % Mu_z (we need to fix this to 1 here as this is
                             % a simple model problem solving
                             % curl curl E - omega^2 E=0
epz = 1;                     % Ep_z
domega =0.01;                 % plane wave component of wave vector
omega = sqrt(domega^2+domega^2+domega^2); % Omega

% Mu_r, Ep_r, Sigma, J
% Mat 1
mu(1) = 1;
epl(1) = 1;
sigma(1) = 0.;
jsrc(1,1:3) = [0  0  0];
% Mat 2
mu(2) = 1;
epl(2) = 1;
sigma(2) = 0;
jsrc(2,1:3) = [0  0  0];

% % for head model
% for i=1:18
% mu(i)=1;
% epl(i)=1;
% jsrc(i,1:3)=[0 0 0];
% sigma(i)=0;
% end
delta = 1;  % Object Size

probdata.matr.nmat=nmat;
probdata.matr.muz=muz;
probdata.matr.epz=epz;
probdata.matr.omega=omega;
probdata.matr.mu=mu;
probdata.matr.epl=epl;
probdata.matr.sigma=sigma;
probdata.matr.jsrc=jsrc;
probdata.matr.delta=delta;
probdata.matr.shift=[0;0;0]';
probdata.matr.matcond=[1 2]; % specify regions where gradients basis 
                             % functions to be included
% probdata.matr.matcond=[1:18];
%-------------------------------------------------------------material data
probdata.sol.regopt=3;       % indefinite system so use direct solve
% define boundary dat

% current bctype
% bctype 5 Object (no bcs)
% bctype 3 Neumann type
% bctype 2 Dirichlet

arg.domega=domega;

% define the function handle for Dirichlet BC's
probdata.es.dirfun=@esproblemdir;
probdata.es.dirfunarg=arg;

% define the source terms
probdata.es.srcfunarg=[];
probdata.es.srcfun=@esproblemsrc;

% % define the function handle for Neumann BC's
probdata.es.neufun=@esproblemneu;
probdata.es.neufunarg=arg;

% exact solution
probdata.es.exactfun=@esproblemexact;
probdata.es.exactfunarg=arg;

probdata.es.exactcurlfun=@esproblemexactcurl;
probdata.es.exactcurlfunarg=arg;

%--------------------------------------------------------------------------

% Dirichlet Boundary Condition---------------------------------------------
function e=esproblemdir(x,y,z,index,arg)       
e=zeros(3,1);
domega = arg.domega;
if index==2
    % this would be non-zero for non zero type BC's
    exp=y+x+z;
    e(1)=complex(cos(domega*exp),sin(domega*exp));
    e(2)=-complex(cos(domega*exp),sin(domega*exp));
    e(3)=complex(0,0);
else
    e(1)=complex(0,0);
    e(2)=complex(0,0);
    e(3)=complex(0,0);
end

% Neumann Boundary Condition-----------------------------------------------
function [curle]=esproblemneu(x,y,z,index,arg)         % Exact curl field
curle=zeros(3,1);
domega = arg.domega;
exp=y+x+z;
curle(1)=domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));
curle(2)=domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));
curle(3)=-2*domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));



function curle=esproblemcurl(x,y,z,arg)
curle=zeros(3,1);
domega = arg.domega;
exp=y+x+z;
curle(1)=domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));
curle(2)=domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));
curle(3)=-2*domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));

% Src term-----------------------------------------------------------------

function src=esproblemsrc(x,y,z,dum,arg,domain)               % Source term for es
src=zeros(3,1);
% Analytical Solution Field------------------------------------------------

function e=esproblemexact(x,y,z,arg)               % Exact field
e=zeros(3,1);
domega = arg.domega;
exp=y+x+z;
e(1)=complex(cos(domega*exp),sin(domega*exp));
e(2)=-complex(cos(domega*exp),sin(domega*exp));
e(3)=complex(0,0);

% Analytical Solution Curl Field-------------------------------------------

function [curle,dum]=esproblemexactcurl(x,y,z,arg)         % Exact curl field
curle=zeros(3,1);
domega = arg.domega;
exp=y+x+z;
curle(1)=domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));
curle(2)=domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));
curle(3)=-2*domega*complex(0,1)*complex(cos(domega*exp),sin(domega*exp));
dum=[];
%--------------------------------------------------------------------------
