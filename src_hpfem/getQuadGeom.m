% Function builds coefficents for geometry approximation using an L2
% projection from quadratic mesh file to a set of heirarchic polynomials of
% degree gorder+1

function [edgecof,facecof,mesh]=getQuadGeom(mesh,gorder,probdata)

% Extract relevant data from mesh structure
nelem=mesh.Nelements;
cintma=mesh.cintma;
intma=mesh.intma;
coord=mesh.Coordinates;
nside=mesh.edge.nside;
nface=mesh.face.nface;

mycoord=zeros(10,3,nelem);
if gorder ~=0
    % loop over elements in mesh
    for i =1:nelem 
                
                % transfer coordinates to local array
                xy = coord(intma(i,1:4),1:3);
                xy(5:10,1:3) = coord(cintma(i,1:6),1:3);
                mycoord(1:10,1:3,i) = xy(1:10,1:3);
                         
    end        
end


mesh.mycoord=mycoord;
mesh.flagBlend=ones(nelem,1);

% zero edge coefficents
edgecof = zeros(nside,3*(gorder-1)+3);

% zero face coefficents
facecof = zeros(nface,(gorder*(gorder-1)/2-1)*3+3);
