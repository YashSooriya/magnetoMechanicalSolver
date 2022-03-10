% Function to choose reference element type for each tetrahedron in the
% mesh, also renumbers vertices to reflect one of two standard types.

function [mesh]=gettype(mesh)

% Extract relevant data from mesh structure
nelem=mesh.Nelements;
coord=mesh.Coordinates;
intma=mesh.intma;
cintma=mesh.cintma;
quadlin=mesh.quadlin;

% Initialize variables
new = zeros(1,4);
oldedge = zeros(6,2);
newedge = zeros(6,2);
cold = zeros(1,6);
v = zeros(4,3);
xy = zeros(10,3);
eltype = zeros(1,nelem);



% plot out an element
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


for j = 1:4
    for k = 1:3
        xy(j,k)=coord(intma(20,j),k);
    end
end

if quadlin==2
for j = 1:6
    for k = 1:3
        xy(j+4,k)=coord(cintma(20,j),k);
    end
end

for i = 1:nelem
    cold = cintma(i,:);
    
    cintma(i,1) = cold(1);
    cintma(i,2) = cold(4);
    cintma(i,3) = cold(2);
    cintma(i,4) = cold(3);
    cintma(i,5) = cold(5);
    cintma(i,6) = cold(6);
end
end

for i =1:nelem
    if quadlin ==2
        oldedge(1,1) = intma(i,1);
        oldedge(1,2) = intma(i,2);
        oldedge(2,1) = intma(i,2);
        oldedge(2,2) = intma(i,3);
        oldedge(3,1) = intma(i,3);
        oldedge(3,2) = intma(i,1);
        oldedge(4,1) = intma(i,1);
        oldedge(4,2) = intma(i,4);
        oldedge(5,1) = intma(i,2);
        oldedge(5,2) = intma(i,4);
        oldedge(6,1) = intma(i,3);
        oldedge(6,2) = intma(i,4);
        cold = cintma(i,:);
    end
    
% set the max and min indicies in an element and store their location        
    maxl = 0;
    minl = 10000000;
    for j = 1:4
        if maxl <intma(i,j)
            maxl = intma(i,j);
            flagmax=j;
        end
        if minl>intma(i,j)
            minl = intma(i,j);
            flagmin = j;
        end
    end
    
    new(1)=minl;
    new(4)=maxl;
    k = 1;
    for j=1:4
        if j ~= flagmin && j ~=flagmax
            k = k+1;
            new(k) = intma(i,j);
        end
    end
    
    for j = 1:4
        xy(j,1:3) = coord(new(j),1:3);
    end
    
% compute v1,v2,v3 and check orinetation of tet    
    v1 = xy(2,:) - xy(1,:);
    v2 = xy(3,:) - xy(1,:);
    v3 = xy(4,:) - xy(1,:);
    
% check scaler tripple product v3.v1 x v2 >0
    tp = scvectp(v3,v1,v2);
    if tp >0
        intma(i,1:4) = new(1:4);
    else
        k = new(2);
        new(2) = new(3);
        new(3) = k;
        intma(i,1:4) = new(1:4);
    end
    if intma(i,2) < intma(i,3) 
        eltype(i) = 1;
    else
        eltype(i) = 2;
    end
    
    if quadlin == 2 
        newedge(1,1) = intma(i,1);
        newedge(1,2) = intma(i,2);
        newedge(2,1) = intma(i,2);
        newedge(2,2) = intma(i,3);
        newedge(3,1) = intma(i,3);
        newedge(3,2) = intma(i,1);
        newedge(4,1) = intma(i,1);
        newedge(4,2) = intma(i,4);
        newedge(5,1) = intma(i,2);
        newedge(5,2) = intma(i,4);
        newedge(6,1) = intma(i,3);
        newedge(6,2) = intma(i,4);
        
        flag = ones(1,6);
        
        for j = 1:6
            for k = 1:6
                if (newedge(j,1)==oldedge(k,1) && newedge(j,2)==oldedge(k,2)) || (newedge(j,2)==oldedge(k,1) && newedge(j,1)==oldedge(k,2))
                    flag(j) = 0;
                    cintma(i,j) = cold(k);
                end
            end
        end
        
        for j = 1:6
            if flag(j) ~= 0
                error(message('stop'));
            end
        end
    end    
end

mesh.intma=intma;
mesh.cintma=cintma;
mesh.eltype=eltype;
display('compute orinentation types for all elements');
