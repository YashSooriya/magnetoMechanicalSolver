clear all
clc

solver = 1;             % Use parallel solver (1) or serial solver (0)
orderEM = 1;            % Set polynomial order of the H(curl) basis functions
orderMech = 1;          % Set polynomial order of the H^1 basis functions
Ns = 90;                % Set number of snapshots desired
axis = 'z';             % Set non-zero dirichlett direction (x,y,z)
dampRatio = 2e-3;          % Set damping ratio
Ncores = 8;             % Set number of cores to assign to problem in parallel solver


disp_even = 0;          % Use even displacement distribution (1) or specified displacements (0)

disp_min = 0;           % Set starting displacement value (in metres)
disp_max = 0.1;         % Set final displacement value (in metres)
disp_inc = 0.02;        % Set increment steps (in metres)

disp_spec = ([0 0.0005 0.001 0.005 0.02 0.08 0.1]);     % Set specified displacement values


readMesh = 0;           % Initialise readMesh as 1

if solver == 1
    if disp_even == 0
        for i = 1:length(disp_spec)
            if axis == 'x'
                if i == 1
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                else
                    readMesh = 0;
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                end
            elseif axis == 'y'
                if i == 1
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                else
                    readMesh = 0;
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                end
            elseif axis == 'z'
                if i == 1
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh)
                else
                    readMesh = 0;
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh)
                end
            else
                disp('Axis for non-0 dir has been defined incorrectly')
            end
        end
    elseif disp_even == 1
        for i = disp_min:disp_inc:disp_max
            if axis == 'x'
                if i == 1
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                else
                    readMesh = 0;
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                end
            elseif axis == 'y'
                if i == 1
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                else
                    readMesh = 0;
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh)
                end
            elseif axis == 'z'
                if i == 1
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh)
                else
                    readMesh = 0;
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh)
                end
            else
                disp('Axis for non-0 dir has been defined incorrectly')
            end
        end
    end
elseif solver == 0
    for i = disp_min:disp_inc:disp_max
        if axis == 'x'
            mainSerialBatch(1,1,Ns,[i,0,0])
        elseif axis == 'y'
            mainSerialBatch(1,1,Ns,[0,i,0])
        elseif axis == 'z'
            mainSerialBatch(1,1,Ns,[0,0,i])
        else
            disp('Axis for non-0 dir has been defined incorrectly')
        end
    end
else
    disp('Solver choice incorrect')
end