% Function assigns a unique number to each edge in the mesh and builds the
% element to edge connectivities

% glob = the element to global edge number connectivity.

function [mesh,lspo]=side_mech(mesh)

% Extract relevant data from mesh structure
nelem=mesh.Nelements;
npoin=mesh.nNodes;
eltype=mesh.eltype;
intma=mesh.intma;

glob = zeros(nelem,6);     % glob = the local to global edge connectivities
maxside = 2*(nelem+npoin-1);
lspo = zeros(maxside,2);
lo=1;
iside = zeros(1,2*maxside); % iside = two node numbers for each face


ilo = [2,3;
       1,3;
       1,2;
       1,4;
       2,4;
       3,4];
   
ilo2 = [3,2;
       1,3;
       1,2;
       1,4;
       2,4;
       3,4];

lpo = zeros(1,npoin);

nside = 0;
for ie = 1:nelem %200
    for js = 1:6 %201
        if eltype(ie) ==1
            i1 = intma(ie,ilo(js,1));
            i2 = intma(ie,ilo(js,2));
        else
            i1 = intma(ie,ilo2(js,1));
            i2 = intma(ie,ilo2(js,2));
        end
        
        j1 = i1;
        j2 = i2;
        l1 = lpo(j1);
        
%         if l1 ~=0 && lspo(l1,1) ==j2
%             glob(ie,js) = nside;
%             continue
%         end
 
%         while l1 ~= 0 
%             if lspo(l1,1)==js
%                 glob(ie,js) = nside;
%                 break
%             end
%             lo = l1;
%             l1 = lspo(lo,2);
%         end
        while 1
        if l1==0
            nside = nside + 1;
            glob(ie,js) = nside;
            if nside > maxside
                error(message('nside > maxside'));
            end
            lspo(nside,1) = j2;
            lspo(nside,2) = 0;
            if lpo(j1) ==0
                lpo(j1) = nside;
            else
                lspo(lo,2) = nside;
            end
            break
        else
            if lspo(l1,1) == j2
                glob(ie,js)=l1;
                break
            end
            lo=l1;
            l1 = lspo(lo,2);
        end
        end  
    end %201
end %200

for ip = 1:npoin
    l1 = lpo(ip);
    while l1 ~=0
        iside(l1) = ip;
        iside(l1+nside) = lspo(l1,1);
        l1 = lspo(l1,2);
    end
end

mesh.edge.nside=nside;










