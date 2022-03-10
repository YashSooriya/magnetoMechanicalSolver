function[mesh]=centerz(mesh)
coord=mesh.Coordinates;
coord(:,3)=coord(:,3)-2.5;
mesh.Coordinates=coord;