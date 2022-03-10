function [mesh]=BlendingFunction(mesh,probdata)
gorder=probdata.jb.gorder;
nside=mesh.edge.nside;
nface=mesh.face.nface;
nside_mech=mesh.mech.nside;
nface_mech=mesh.mech.nface;
% zero edge coefficents
edgecof = zeros(nside,3*(gorder-1)+3);
edgecof_mech = zeros(nside_mech,3*(gorder-1)+3);
% zero face coefficents
facecof = zeros(nface,(gorder*(gorder-1)/2-1)*3+3);
facecof_mech = zeros(nface_mech,(gorder*(gorder-1)/2-1)*3+3);

mesh.mech.face.cond=mesh.mech.cond;
mesh.mech.face.nface=mesh.mech.nface;
mesh.mech.face.globfa=mesh.mech.globfa;
mesh.mech.edge.nside=mesh.mech.nside;
mesh.mech.edge.bcedge=mesh.mech.bcedge;
mesh.mech.edge.glob=mesh.mech.glob;



if gorder >=1
    % The currently available curved geometry types
    % 1 = not available
    % 2 = sphere exact geometry
    % 3 = not available
    % 4 = quadratic geometry file
    display(['The geometry type is ',num2str(gorder)]);
    if gorder ==4
        [edgecof_mech,facecof_mech,mesh.mech]=getblendh1(mesh.mech,gorder,probdata);
         [edgecof,facecof,mesh]=getblendh1(mesh,gorder,probdata);

    elseif gorder==2
 
             [edgecof,facecof,mesh]=getblendqh1(mesh,gorder,probdata);
            [edgecof_mech,facecof_mech,mesh.mech]=getblendqh1(mesh.mech,gorder,probdata);
            
    display('completed all bending');
    end
end

% Store blending fuction information in mesh structure
mesh.edgecof=edgecof;
mesh.facecof=facecof;
mesh.mech.edgecof=edgecof_mech;
mesh.mech.facecof=facecof_mech;
