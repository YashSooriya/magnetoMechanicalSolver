clear all
clc

solver = 1;             % Use parallel solver (1) or serial solver (0)
Ns = 90;                % Set number of snapshots desired

axis = 'z';             % Set non-zero dirichlett direction (x,y,z)
<<<<<<< Updated upstream
disp_min = 0;           % Set starting displacement value (in metres)
disp_max = 0.1;         % Set final displacement value (in metres)
disp_inc = 0.02;        % Set increment steps (in metres)
=======
disp_min = 0.0005;           % Set starting displacement value (in metres)
disp_max = 0.005;        % Set final displacement value (in metres)
disp_inc = 0.0045;        % Set increment steps (in metres)
>>>>>>> Stashed changes

if solver == 1
    for i = disp_min:disp_inc:disp_max
        if axis == 'x'
            mainParallelBatch(Ns,[i,0,0])
        elseif axis == 'y'
            mainParallelBatch(Ns,[0,i,0])
        elseif axis == 'z'
            mainParallelBatch(Ns,[0,0,i])
        else
            disp('Axis for non-0 dir has been defined incorrectly')
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
