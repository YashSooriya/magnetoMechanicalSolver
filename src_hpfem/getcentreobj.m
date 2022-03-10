% Function finds the centre of the mesh and moves the coordinates such that
% the object is centered at (0,0,0) 
function [mesh]=getcentreobj(mesh)

% find object centre
mat=mesh.mat;
coord=mesh.Coordinates;
matc=mesh.matc;
type=mesh.eltype;
nelem=mesh.Nelements;
intma=mesh.intma;
npoin=mesh.nNodes;
[cen]=centreobj(coord,mat,matc,type,nelem,intma);

% adjust coordinates
for i=1:npoin
    for j=1:3
        coord(i,j)=coord(i,j)-cen(j);
    end
end

disp('Correct coordinates so that origin is centre of object')
% find object centre
[cen]=centreobj(coord,mat,matc,type,nelem,intma);

mesh.Coordinates=coord;


function [cen]=centreobj(coord,mat,matc,type,nelem,intma)

myflag=0;

vol=0.d0;
cen=[0;0;0];
gorder=0;
gesizet = (gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
maxnip=100;
nipd=ceil(max((2*(1+1)+1)/2,2));
% obtain integration points......over a tetrahedra
[intxi,inteta,intzeta,intw,nip]=intpoints(nipd,maxnip);
% set up mapping functions
axi(1)=1;
axi(2)=0;
axi(3)=0;

aeta(1)=0;
aeta(2)=1;
aeta(3)=0;

azeta(1)=0;
azeta(2)=0;
azeta(3)=1;

asxi(1)=1;
asxi(2)=0;
asxi(3)=0;

aseta(1)=0;
aseta(2)=1;
aseta(3)=0;

aszeta(1)=0;
aszeta(2)=0;
aszeta(3)=1;
%store geometry basis functions over volume (since gorder =0 we don't need
%to worry about type 1/type 2 reference elements)
gphx1=zeros(nip,gesizet);
gphy1=zeros(nip,gesizet);
gphz1=zeros(nip,gesizet);

phh11=zeros(nip,gesizet);

for i=1:nip
    gph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,gorder+1,1,gesizet);
    
    gphx1(i,1:gesizet)=gph(1:gesizet,1)';
    gphy1(i,1:gesizet)=gph(1:gesizet,2)';
    gphz1(i,1:gesizet)=gph(1:gesizet,3)';
    
    phh1=basish1(gesizet,intxi(i),inteta(i),intzeta(i),gorder+1,1);
    phh11(i,1:gesizet)=phh1(1:gesizet,1)';
    
end

% find the volume of the conductive region by integration
% locate also the centre of the conductive region
for i=1:nelem
    if matc(i)==1
        
        xy = coord(intma(i,1:4),1:3);
        
        %transfer coefficents to locations vertices
        mycoord(1:4,1:3) = xy(1:4,1:3);
        
        flag=0;
        le=[];lf=[];
        for pp=1:nip
            
            gph(1:gesizet,1)=gphx1(pp,1:gesizet)';
            gph(1:gesizet,2)=gphy1(pp,1:gesizet)';
            gph(1:gesizet,3)=gphz1(pp,1:gesizet)';
            ph1(1:gesizet,1)=phh11(pp,1:gesizet)';
            %             [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian(xy,intxi(pp),inteta(pp),intzeta(pp),...
            %                 type,le,lf,flag,gorder,mycoord);
            [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord) ;
            
            %             [x,y,z]=getxyz(xy,intxi(pp),inteta(pp),intzeta(pp));
            [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
            
            cen(1) =cen(1) +x *det *intw(pp);
            cen(2) =cen(2) +y *det *intw(pp);
            cen(3) =cen(3) +z *det *intw(pp);
            
            vol=vol+det*intw(pp);
            
        end
    end
end

%     normalise coordinates
cen(1)=cen(1)/vol;
cen(2)=cen(2)/vol;
cen(3)=cen(3)/vol;


disp('The centre of the object is')
disp(['(',num2str(cen(1)),',',num2str(cen(2)),',',num2str(cen(3)),')'])

