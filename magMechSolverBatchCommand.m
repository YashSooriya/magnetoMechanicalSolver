clear all
clc

iter = 3;

%-------------------------------------------------------------------------
% Solver Definition
%-------------------------------------------------------------------------
solver = [1; 1; 1];                        % Use parallel solver (1) or serial solver (0)
orderEM = [2; 3; 3];                       % Set polynomial order of the H(curl) basis functions
orderMech = [2; 3; 3];                     % Set polynomial order of the H^1 basis functions
Ns = [90; 90; 90];                         % Set number of snapshots desired
Ncores = [8; 8; 8];                        % Set number of cores to assign to problem in parallel solver
singleStatic = [1; 1; 1];                  % Set if static solver is only used once (1) or not (0)
readMesh = [0; 0; 0];                      % Initialise readMesh as 1

%-------------------------------------------------------------------------
% Problem Definition
%-------------------------------------------------------------------------
axis = ['z'; 'z'; 'z'];                    % Set non-zero dirichlett direction (x,y,z)
dampRatio = [0; 0; 0];                     % Set damping ratio

%-------------------------------------------------------------------------
% Displacement Definition
%-------------------------------------------------------------------------
disp_even = [0; 0; 0];                     % Use even displacement distribution (1) or specified displacements (0)

disp_min = [0; 0; 0];                      % Set starting displacement value (in metres)
disp_max = [0.1; 0.1; 0.1];                % Set final displacement value (in metres)
disp_inc = [0.1; 0.1; 0.1];                % Set increment steps (in metres)

dir_DC = [0; 0; 0];                        % Set if Dirichlett displacement should be non-zero (1) or not (0) for DC solver

%-------------------------------------------------------------------------
% Output Frequency Refinement Definition
%-------------------------------------------------------------------------
refined = [0; 0; 0];                       % Set if frequency out increment is refined (1) or not (0)
 
POIs = [1000; 1000; 1000];                 % Set points of interest for refinement
range = [100; 100; 100];                   % Set range of frequencies before and after POIs to be refined
refinedDelFOut = [5; 5; 5];                % Set refined frequency increment

%-------------------------------------------------------------------------
% Postprocessing
%-------------------------------------------------------------------------
method = [1; 0; 1];                        % Use updated faster power and energy computation method (1) or not (0)


for i = 1:iter
    if i == 10
        disp_spec = ([0]);
    elseif i == 20
        disp_spec = ([0]);
    elseif i == 30
        disp_spec = ([0.001 0.005 0.02 0.08 0.1]);
    elseif i == 1
        disp_spec = ([0.001 0.005]);
    elseif i == 2
        disp_spec = ([0]);
    elseif i == 3
        disp_spec = ([0.001 0.005 0.02 0.08 0.1]);
    end
    magMechSolverBatch(solver(i), orderEM(i), orderMech(i), Ns(i), Ncores(i), singleStatic(i), ...
    readMesh(i), axis(i), dampRatio(i), disp_even(i), disp_min(i), disp_max(i), disp_inc(i), disp_spec, ...
    dir_DC(i), refined(i), POIs(i), range(i), refinedDelFOut(i), method(i))
end