function [coordn,intman,Splitref,splitt]=split(nelem,npoin,nside,nface,coord,intma,glob,globfa,order,edgecof,facecof,gorder,eltype)

%[coordn,intman,splitt,v]=split(nelem,npoin,nside,nface,coord,intma,glob,globfa)

% % ped = no splitting
% % ped =1 one point per edge
% % ped =2 two points per edge with one ponts per face
% ped = 0;
% 
% % reference element coord
% v(1,1)=-1;
% v(1,2)=0;
% v(1,3)=0;
% 
% v(2,1)=1;
% v(2,2)=0;
% v(2,3)=0;
% 
% v(3,1)=0;
% v(3,2)=sqrt(3);
% v(3,3)=0;
% 
% v(4,1)=0;
% v(4,2)=sqrt(3)/3;
% v(4,3)=2*(sqrt(2)/sqrt(3));

% produce a splitting of the reference element
Splitref = nodalCoordinatesTetraMinus1(order+2);
[nptet dum]=size(Splitref);

% generate a Delaunay triangulation of this tetrahedron
intmasplit = delaunay(Splitref(:,1),Splitref(:,2),Splitref(:,3));
[netet dum]=size(intmasplit);

% produce a mesh for plotting consisting of split tetrahedrons
npoinn=0;
intman=[];
lec=[];
lfc=[];
for i=1:nelem
    % generate connectivities
    intman=[intman; npoinn+intmasplit];
    
    xy = coord(intma(i,1:4),1:3);
    
    flag=0;
    for j=1:6
        for p=1:gorder
            for k=1:3
                lec(j,p,k)=edgecof(glob(i,j),((p-1)*3)+k);
                if abs(edgecof(glob(i,j),((p-1)*3)+k))>0.0000001
                    flag=1;
                end
            end
        end
    end
    
    for j=1:4
        for p=1:gorder*(gorder-1)/2
            for k=1:3
                lfc(j,p,k)=facecof(globfa(i,j),((p-1)*3)+k);
            end
        end
    end
    
    % generate coordinates
    for j=1:nptet
        % computx x,y,z
        [x,y,z]= getxyzcu(xy,Splitref(j,1),Splitref(j,2),Splitref(j,3),lec,lfc,flag,gorder,eltype(i));
        %[x,y,z]=getxyz(xy,Splitref(j,1),Splitref(j,2),Splitref(j,3));
        npoinn=npoinn+1;
        coordn(npoinn,:) = [x,y,z];
        splitt(i,j)=npoinn;
    end
end
    
    


% % define number of points per element depending on the splitting
% if ped==0
%     nppe=4;
% elseif ped ==1
% 
% nppe = 4+6;
% coordn = zeros(npoin+nside,3);
% splitt = zeros(nelem,nppe);
% 
% v(5,:) = (v(3,:)+v(2,:))/2; % edge 1
% v(6,:) = (v(3,:)+v(1,:))/2; % edge 2
% v(7,:) = (v(2,:)+v(1,:))/2; % edge 3 
% v(8,:) = (v(4,:)+v(1,:))/2; % edge 4 
% v(9,:) = (v(4,:)+v(2,:))/2; % edge 5 
% v(10,:) = (v(4,:)+v(3,:))/2; % edge 6 
% 
% elseif ped==2
% 
% nppe = 4+2*6+4;
% coordn = zeros(npoin+2*nside+nface,3);
% splitt = zeros(nelem,nppe);
% 
% % edge 1
% v(5,:) =  (v(3,:)-v(2,:))/3+v(2,:); 
% v(6,:) =  2*(v(3,:)-v(2,:))/3+v(2,:);
% 
% % edge 2
% v(7,:) = (v(3,:)-v(1,:))/3+v(1,:); 
% v(8,:) = 2*(v(3,:)-v(1,:))/3+v(1,:);
% % edge 3 
% v(9,:) = (v(2,:)-v(1,:))/3+v(1,:); 
% v(10,:) = 2*(v(2,:)-v(1,:))/3+v(1,:);
% % edge 4 
% v(11,:) = (v(4,:)-v(1,:))/3+v(1,:); 
% v(12,:) = 2*(v(4,:)-v(1,:))/3+v(1,:); 
% % edge 5 
% v(13,:) = (v(4,:)-v(2,:))/3+v(2,:); 
% v(14,:) = 2*(v(4,:)-v(2,:))/3+v(2,:); 
% % edge 6 
% v(15,:) = (v(4,:)-v(3,:))/3+v(3,:); 
% v(16,:) = 2*(v(4,:)-v(3,:))/3+v(3,:);
% 
% % face 1
% v(17,:) = (v(2,:)+v(3,:)+v(4,:))/3;
% % face 2
% v(18,:) = (v(1,:)+v(3,:)+v(4,:))/3;
% % face 3
% v(19,:) = (v(1,:)+v(2,:)+v(4,:))/3;
% % face 3
% v(20,:) = (v(1,:)+v(2,:)+v(3,:))/3;
% 
% end
% 
% coordn(1:npoin,:) = coord;
% 
% for i = 1:nelem
%     xy = coord(intma(i,1:4),1:3);
%     
%     splitt(i,1:4) = intma(i,1:4);
%     
%     if ped ==1
%         for j = 1:6 % for each edge
%             xi = v(4+j,1);
%             eta = v(4+j,2);
%             zeta = v(4+j,3);
%             [x,y,z]=getxyz(xy,xi,eta,zeta);
%             coordn(npoin+ glob(i,j),:) = [x,y,z];
%             
%             
%             
%             splitt(i,4+j) = npoin+ glob(i,j);
%         end
%     elseif ped==2
%         for j = 1:6 % for each edge
%             xi = v(4+2*j-1,1);
%             eta = v(4+2*j-1,2);
%             zeta = v(4+2*j-1,3);
%             [x,y,z]=getxyz(xy,xi,eta,zeta);
%             coordn(npoin+ glob(i,j),:) = [x,y,z];
%             splitt(i,4+2*j-1) = npoin+ 2*glob(i,j)-1;
%             
%             xi = v(4+2*j,1);
%             eta = v(4+2*j,2);
%             zeta = v(4+2*j,3);
%             [x,y,z]=getxyz(xy,xi,eta,zeta);
%             coordn(npoin+ nside+glob(i,j),:) = [x,y,z];
%             splitt(i,4+2*j) = npoin+ 2*glob(i,j);
%         end
%         for k = 1:4 % for each face
%             xi = v(16+k,1);
%             eta = v(16+k,2);
%             zeta = v(16+k,3);
%             [x,y,z]=getxyz(xy,xi,eta,zeta);
%             coordn(npoin+ 2*nside+globfa(i,k),:) = [x,y,z];
%             
%             splitt(i,16+k) = npoin+ 2*nside+globfa(i,k);
%         end
% 
%     end
% end
% if ped==1 || ped==2
%     intman = delaunay(coordn(:,1),coordn(:,2),coordn(:,3));
% else 
%     intman=intma;
% end