% Set up Standard Options for BC's, src data etc for an MPT calculation
function probdata=standoptions(probdata)


%-------------------------------------------------------------------------
% Not used
% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 0; 

% % plot the  H_alpha-H_0 - D2G PT H_0  on a line
% % darkooption = 0 - do not plot
% % darkooption = 1 - plot
 darkooption = 0; 

% % plot the  voltage  on a line
% % voltoption = 0 - do not plot
% % voltoption = 1 - plot
 voltoption = 0; 
 
probdata.mesh.voltoption=voltoption;
probdata.mesh.darkooption=darkooption;
 


% solver options ---------------------------------------------------------

probdata.sol.regopt =2 % 1- use regularisation with cg for gradent blocks
                       % 2- use regularisation with direct solve for gradient blocks
                       % 3- direct solve

                       
probdata.jb.nrhs=1;  % 3 right hand sides
%--------------------------------------------------------------------------
% Standard mesh options
matflg = 1;        % Mat Flag from file
ppus   = 1 ;       % Parrellel Processors (not used)
% Polynomial Degree Info
ptype = 2;         % (Type)

probdata.jb.matflg=matflg;
probdata.jb.ppus = ppus;
probdata.jb.ptype = ptype; 

%-------------------------------------------------------------------------
% Standard Geometry Options

% Blending Function Info
gorder =2;        % order (if gorder > 0, quadlin =2 required)
g1 = 4;            % type (g1=4 recomended)
sufv(1)= 5;        % surfaceval
gag = 2;           % Goagain  % 1- go again 2- once
svchk = 0;         % Do check surface area /volume as not expecting a sphere

probdata.jb.gorder = gorder;
probdata.jb.g1 = g1;
probdata.jb.sufv = sufv(1);
probdata.jb.gag = gag(1);
probdata.mesh.svchk=svchk;
%-------------------------------------------------------------------------

