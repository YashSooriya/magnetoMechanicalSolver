function [Static,UnknownStatic]=staticSolver(Mesh,Basis,Quadrature,ProblemData,Options)

disp('--------------------------------------------')
disp('Running the static solver')
disp('--------------------------------------------')

%=========================================================================
% Compute the static field
%=========================================================================
% Initialise the frequency (rad/s) - 1 due to regularisation
omega=1;
ProblemData.matr.omega=omega;

% Initialise the static analysis flag
ProblemData.probstatic=1;

% Define the K, C and M weights 
w0                     = 1;
w1                     = 1;
w2                     = 0;

[UnknownStatic]=elementUnknownNumberingStatic(Mesh,ProblemData);
%-------------------------------------------------------------------------
% Initial guesses of the electromagnetic field
%-------------------------------------------------------------------------

% Define the problem flag
ProblemData.probFlag=1;

% Initial guess of the solution field
[ADC] = initialGuess(Mesh,Basis,Quadrature,UnknownStatic.EM,ProblemData,omega);

%-------------------------------------------------------------------------
% Initial guess of the mechanical field
%--------------------------------------------------------------------------

% Define the problem flag
ProblemData.probFlag=2;

% Initial guess for the mechanical problem
[UDC] = initialGuess(Mesh,Basis,Quadrature,UnknownStatic.Mech,ProblemData,omega);


%-------------------------------------------------------------------------
% Store the solutions
%-------------------------------------------------------------------------
% Compute the initial solution vector
X                      = [ADC;UDC];
dXdt                   = w1*X;
dXdt2                  = w2*X;

% Store the solutions of the static field in the static data structure
Static.sol             = [X dXdt dXdt2];
%-------------------------------------------------------------------------
% Monolithic solver
%-------------------------------------------------------------------------
% Run the monolithic solver

StaticCurrent=Static;
%U=zeros(length(X),1);
% Solve the electromagnetic problem

[U,~,UnknownStatic] = linearSystemSolver1(Static,StaticCurrent,Mesh,Basis,Quadrature,UnknownStatic,ProblemData,[X dXdt dXdt2],w0,w1,w2,Options);


%--------------------------------------------------------------------------
% Update solution (converges in one iteration)
X=X+U;
Static.sol=[X dXdt dXdt2];
% Solve mechanical static problem


%-------------------------------------------------------------------------
% Store the solution field data
%-------------------------------------------------------------------------
% Store the results in the statics structure
Static.sol=[X dXdt dXdt2];
%Static.conditionNum = condNum;
Static.nSparse=UnknownStatic.system.nSparse;


end