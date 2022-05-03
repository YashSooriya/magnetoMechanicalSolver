clear all
clc

%-------------------------------------------------------------------------
% Solver Definition
%-------------------------------------------------------------------------
solver = 1;                 % Use parallel solver (1) or serial solver (0)
orderEM = 1;                % Set polynomial order of the H(curl) basis functions
orderMech = 1;              % Set polynomial order of the H^1 basis functions
Ns = 90;                    % Set number of snapshots desired
Ncores = 8;                 % Set number of cores to assign to problem in parallel solver
singleStatic = 1;           % Set if static solver is only used once (1) or not (0)
readMesh = 0;               % Initialise readMesh as 1

%-------------------------------------------------------------------------
% Problem Definition
%-------------------------------------------------------------------------
axis = 'z';                 % Set non-zero dirichlett direction (x,y,z)
dampRatio = 0;              % Set damping ratio

%-------------------------------------------------------------------------
% Displacement Definition
%-------------------------------------------------------------------------
disp_even = 0;              % Use even displacement distribution (1) or specified displacements (0)

disp_min = 0;               % Set starting displacement value (in metres)
disp_max = 0.1;             % Set final displacement value (in metres)
disp_inc = 0.02;            % Set increment steps (in metres)

disp_spec = ([0 0.001]);    % Set specified displacement values

dir_DC = 1;                 % Set if Dirichlett displacement should be non-zero (1) or not (0) for DC solver

%-------------------------------------------------------------------------
% Output Frequency Refinement Definition
%-------------------------------------------------------------------------
refined = 0;

POIs = [1000 2000 3500];
range = 100;
refinedDelFOut = 5;

%-------------------------------------------------------------------------
% Start Batch Processing
%-------------------------------------------------------------------------
if solver == 1
    if disp_even == 0
        for i = 1:length(disp_spec)
            if axis == 'x'
                if singleStatic == 1
                    if i == 1
                        mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut)
                    else
                        postStaticSolver([disp_spec(i), 0, 0], dampRatio, orderEM, orderMech)
                    end
                else
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut)
                end
            elseif axis == 'y'
                if i == 1
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut)
                else
                    postStaticSolver([0, disp_spec(i), 0], dampRatio, orderEM, orderMech)
                end
            elseif axis == 'z'
                if i == 1
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut)
                else
                    postStaticSolver([0, 0, disp_spec(i)], dampRatio, orderEM, orderMech)
                end
            else
                disp('Axis for non-0 dir has been defined incorrectly')
            end
        end
    elseif disp_even == 1
        for i = disp_min:disp_inc:disp_max
            if axis == 'x'
                if i == 1
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut)
                else
                    postStaticSolver([disp_spec(i), 0, 0], dampRatio, orderEM, orderMech)
                end
            elseif axis == 'y'
                if i == 1
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut)
                else
                    postStaticSolver([0, disp_spec(i), 0], dampRatio, orderEM, orderMech)
                end
            elseif axis == 'z'
                if i == 1
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut)
                else
                    postStaticSolver([0, 0, disp_spec(i)], dampRatio, orderEM, orderMech)
                end
            else
                disp('Axis for non-0 dir has been defined incorrectly')
            end
        end
    end
elseif solver == 0
    if disp_even == 0
        for i = 1:length(disp_spec)
            if axis == 'x'
                if i == 1
                    mainSerialBatch(1, 1, Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                else
                    readMesh = 0;
                    mainSerialBatch(1, 1, Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                end
            elseif axis == 'y'
                if i == 1
                    mainSerialBatch(1, 1, Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                else
                    readMesh = 0;
                    mainSerialBatch(1, 1, Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                end
            elseif axis == 'z'
                if i == 1
                    mainSerialBatch(1, 1, Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                else
                    readMesh = 0;
                    mainSerialBatch(1, 1, Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                end
            else
                disp('Axis for non-0 dir has been defined incorrectly')
            end
        end
    elseif disp_even == 1
        for i = disp_min:disp_inc:disp_max
            if axis == 'x'
                if i == 1
                    mainSerialBatch(1, 1, Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                else
                    readMesh = 0;
                    mainSerialBatch(1, 1, Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                end
            elseif axis == 'y'
                if i == 1
                    mainSerialBatch(1, 1, Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                else
                    readMesh = 0;
                    mainSerialBatch(1, 1, Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                end
            elseif axis == 'z'
                if i == 1
                    mainSerialBatch(1, 1, Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                else
                    readMesh = 0;
                    mainSerialBatch(1, 1, Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, readMesh, dir_DC)
                end
            else
                disp('Axis for non-0 dir has been defined incorrectly')
            end
        end
    end
else
    disp('Solver choice incorrect')
end