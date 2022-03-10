function [nelb,npb,coordf,out,help1]=pelem(npoin,nface,nelem,coord,intma,...
    cond,globfa,orderel)

no = [4, 2, 3, 0;
      4, 3, 1, 0;
      4, 1, 2, 0;
      1, 3, 2, 0];
% no2 = [3, 4, 2, 0;
%        1, 3, 4, 0;
%        1, 2, 4, 0;
%        1, 2, 3, 0];

v = zeros(4,3);
out = zeros(nface*3,5);
coordf = zeros(nface*3,3);

v(1,1) = -1;
v(1,2) = 0;
v(1,3) = 0;

v(2,1) = 1;
v(2,2) = 0;
v(2,3) = 0;

v(3,1) = 0;
v(3,2) = sqrt(3);
v(3,3) = 0;

v(4,1) = 0;
v(4,2) = sqrt(3)/3;
v(4,3) = 2*(sqrt(2)/sqrt(3));

help2 = zeros(nface,1);
help1 = zeros(nface,3);
help3 = zeros(npoin,1);
nelb = 0;
npb = 0;

% determine which points are on the boundary
for i=1:nelem
    %transter coordinates
    xy = coord(intma(i,1:4),1:3);
    
    %determine node numbers and coordinates
    for j = 1:4
        if cond(i,j)==2
            % increase element number
            nelb = nelb+1;
            help2(globfa(i,j)) = nelb;
            
            % find output value
            orderi = orderel(i);
            
            % determine node numbers inside this element
            for k = 1:3
                npb = npb+1;
                help3(intma(i,no(j,k))) = npb;
                help1(nelb,k) = npb;
                
                % determine coordinates, no(j,k) contains the kth vertex on face j
                % compute x,y,z
                [x,y,z]=getxyz(xy,v(no(j,k),1)*0.999,v(no(j,k),2)*0.999,v(no(j,k),3)*0.999);
                
                coordf(npb,1) = x;
                coordf(npb,2) = y;
                coordf(npb,3) = z;
                
                out(npb,1:5) = orderi*ones(1,5);
            end
        end
    end   
end
display('finished plot of non-uniform p...compute nzeros')