function [Splitref,intmasplit] = nodalCoordinatesTetraMinus1(nDeg)
%
% Xi = nodalCoordinatesTetraMinus1(nDeg)
%
% Define an equally spaced nodal distribution in the reference tetrahedron
% [-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1]
%

h0 = 2/nDeg;
% Four vertexs
Xi = [-1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1];
% First edge
for i=1:nDeg-1
    Xi = [Xi; -1+i*h0  -1 -1];
end
% Face with zeta=-1
for j=1:nDeg-1
    for i=0:nDeg-j
        Xi = [Xi;-1+i*h0  -1+j*h0 -1];
    end
end
% Levels increasing zeta
for k=1:nDeg-1
    for j=0:nDeg-k
        for i=0:nDeg-k-j
            Xi = [Xi;-1+i*h0  -1+j*h0 -1+k*h0];
        end
    end
end
[nptet dum]=size(Xi);
% Generate a Delaunay triangulation of this tetrahedron
% Note that we prefer to do a splitting of this Tet as it will have
% nicer points which MATLAB's delaunayn can handle, if we use
% new points generated below delaunayn can fail to generate all
% required tets and there may be gaps in the mesh!
intmasplit = delaunayn(Xi);
[netet dum]=size(intmasplit);
pen=zeros(nptet,9);
phen=pen;
% vtuk_puvw_write ( 'rubentet.vtu', nptet, netet, ...
%    Xi, intmasplit, pen, phen )
%  WriteToVTK(Xi,intmasplit,phen,pen,'rubentet.vtk')


% Map these Coordinates to my reference Tetrahedron!
[npoints dum]=size(Xi);

% coordinates of my elements
coord=[-1 0 0;
       1 0 0;
       0 sqrt(3) 0;
       0 sqrt(3)/3 2*sqrt(2)/sqrt(3)];
   
% adjust nodal coordinates (so points are just inside the reference tetrahedron to avoid any issues with
% basis functions which might involve division by 0 at the vertices
centroid=zeros(1,3);
for i=1:4
centroid=centroid+(coord(i,:));
end
centroid=1/4*centroid;
coord=0.999*coord+0.001*ones(4,1)*centroid;



Splitref=zeros(npoints,3);
for i=1:npoints
  x=Xi(i,1);
  y=Xi(i,2);
  z=Xi(i,3);
  % Define area coordinates
  n(1)=1/2*(-1-x-y-z);
  n(2)=1/2*(1+x);
  n(3)=1/2*(1+y);
  n(4)=1/2*(1+z);
  
  for j=1:4
  Splitref(i,:)=Splitref(i,:)+n(j)*coord(j,:);
  end
end

pen=zeros(nptet,9);
phen=pen;

% vtuk_puvw_write ( 'paultet.vtu', nptet, netet, ...
%    Splitref, intmasplit, pen, phen )
% WriteToVTK(Splitref,intmasplit,phen,pen,'paultet.vtk')
