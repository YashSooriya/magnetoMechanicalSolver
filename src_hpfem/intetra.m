function [in,onf]= intetra(V1,V2,V3,V4,P)

% Check a point in or ouf of tethrahedron
%
% input: V1,V2,V3,V4  Vertex of tethrahedron Vi = (xi, yi, zi);
%        P            Point need to be checked P = (x, y, z)
% Output: in          1  inside
%                     0  on the face
%                     -1 outside
%          onf        the face of point lies

% refer to:  http://steve.hollasch.net/cgindex/geometry/ptintet.html
% V1 = [1,0,0];
% V2 = [0,1,0];
% V3 = [0,0,0];
% V4 = [0,0,1];
% P = [0.2,0.2,0.2];

d0 = [V1 1;
      V2 1;
      V3 1;
      V4 1];

d1 = [P  1;
      V2 1;
      V3 1;
      V4 1];
 
d2 = [V1 1;
      P  1;
      V3 1;
      V4 1]; 

d3 = [V1 1;
      V2 1;
      P  1;
      V4 1];
  
d4 = [V1 1;
      V2 1;
      V3 1;
      P  1];
  
D0 = det(d0);
D1 = det(d1);
D2 = det(d2);
D3 = det(d3);
D4 = det(d4);
D = [D1,D2,D3,D4];
if D0 ==0
    error('its not tetra');
elseif any(D ==0)
    if D(1)==0 && D(2)==0 && D(3) ==0 && sign(D4)==sign(D0) ||...
       D(1)==0 && D(3)==0 && D(4) ==0 && sign(D2)==sign(D0) ||...
       D(1)==0 && D(2)==0 && D(4) ==0 && sign(D3)==sign(D0) ||...
       D(2)==0 && D(3)==0 && D(4) ==0 && sign(D1)==sign(D0)
        in = 0;
        onf = 0;   
    elseif D(1)==0 && D(2)==0 && sign(D3)==sign(D0) && sign(D4)==sign(D0) ||...
           D(1)==0 && D(3)==0 && sign(D2)==sign(D0) && sign(D4)==sign(D0) ||...
           D(1)==0 && D(4)==0 && sign(D2)==sign(D0) && sign(D3)==sign(D0) ||...
           D(2)==0 && D(3)==0 && sign(D1)==sign(D0) && sign(D4)==sign(D0) ||...
           D(2)==0 && D(4)==0 && sign(D1)==sign(D0) && sign(D3)==sign(D0) ||...
           D(3)==0 && D(4)==0 && sign(D1)==sign(D0) && sign(D2)==sign(D0)
        in = 0;
        onf = 0; 
    elseif D(1) ==0 && sign(D2)==sign(D0) && sign(D3)==sign(D0) && sign(D4)==sign(D0)
        in = 0;
        onf = 1;        
    elseif D(2)==0 && sign(D1)==sign(D0) && sign(D3)==sign(D0) && sign(D4)==sign(D0)
        in = 0;
        onf = 2;
    elseif D(3)==0 && sign(D1)==sign(D0) && sign(D2)==sign(D0) && sign(D4)==sign(D0)
        in = 0;
        onf = 3;
    elseif D(4)==0 && sign(D1)==sign(D0) && sign(D2)==sign(D0) && sign(D3)==sign(D0)
        in = 0;
        onf = 4;
    else
        in = -1;
        onf = 0;
    end 
elseif sign(D1)==sign(D0) && sign(D2)==sign(D0) && sign(D3)==sign(D0)...
        && sign(D4)==sign(D0) %&& (D1+D2+D3+D4)==D0
    in = 1;
    onf = 0;
else
    in = -1;
    onf=0;
end