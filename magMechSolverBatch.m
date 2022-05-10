function magMechSolverBatch(solver,orderEM,orderMech,Ns,Ncores,singleStatic, ...
    readMesh,axis,dampRatio,disp_even,disp_min,disp_max,disp_inc,disp_spec, ...
    dir_DC,refined,POIs,range,refinedDelFOut,method)

%-------------------------------------------------------------------------
% Start Batch Processing
%-------------------------------------------------------------------------
if solver == 1
    if disp_even == 0
        for i = 1:length(disp_spec)
            if axis == 'x'
                if singleStatic == 1
                    if i == 1
                        mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut,method)
                    else
                        postStaticSolver([disp_spec(i), 0, 0], dampRatio, orderEM, orderMech,method)
                    end
                else
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut,method)
                end
            elseif axis == 'y'
                if i == 1
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut,method)
                else
                    postStaticSolver([0, disp_spec(i), 0], dampRatio, orderEM, orderMech,method)
                end
            elseif axis == 'z'
                if i == 1
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut,method)
                else
                    postStaticSolver([0, 0, disp_spec(i)], dampRatio, orderEM, orderMech,method)
                end
            else
                disp('Axis for non-0 dir has been defined incorrectly')
            end
        end
    elseif disp_even == 1
        for i = disp_min:disp_inc:disp_max
            if axis == 'x'
                if i == 1
                    mainParallelBatch(Ns,[disp_spec(i), 0, 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut,method)
                else
                    postStaticSolver([disp_spec(i), 0, 0], dampRatio, orderEM, orderMech,method)
                end
            elseif axis == 'y'
                if i == 1
                    mainParallelBatch(Ns,[0, disp_spec(i), 0], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut,method)
                else
                    postStaticSolver([0, disp_spec(i), 0], dampRatio, orderEM, orderMech,method)
                end
            elseif axis == 'z'
                if i == 1
                    mainParallelBatch(Ns,[0, 0, disp_spec(i)], dampRatio, orderEM, orderMech, Ncores, readMesh, dir_DC, refined, POIs, range, refinedDelFOut,method)
                else
                    postStaticSolver([0, 0, disp_spec(i)], dampRatio, orderEM, orderMech,method)
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