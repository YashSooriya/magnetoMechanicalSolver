clear
clc

axis = 'z';             % Set non-zero dirichlett direction (x,y,z)
disp_min = 0;           % Set starting displacement value (in metres)
disp_max = 0.1;        % Set final displacement value (in metres)
disp_inc = 0.02;        % Set increment steps (in metres)

for i = disp_min:disp_inc:disp_max
    if axis == 'x'
        mainSerialBatch(1,1,90,[i,0,0])
    elseif axis == 'y'
        mainSerialBatch(1,1,90,[0,i,0])
    elseif axis == 'z'
        mainSerialBatch(1,1,90,[0,0,i])
    else
        disp('Axis for non-0 dir has been defined incorrectly')
    end
end
