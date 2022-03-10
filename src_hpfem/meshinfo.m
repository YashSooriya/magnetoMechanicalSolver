% Function reads in the mesh information from a NETGEN .vol file

function [mesh]=meshinfo(meshdata,ProblemData)

% Extract data from ProblemData structure
meshtype=ProblemData.jb.meshtype;
matflg=ProblemData.jb.matflg;
shift=ProblemData.matr.shift;
delta=ProblemData.matr.delta;


nboun = str2num(meshdata{9});
bsido = zeros(nboun,6);
nelem = str2num(meshdata{9+nboun+5});
help1 = zeros(2*nboun,5);
intma = zeros(nelem,4);
cintma = zeros(nelem,6);
mat = zeros(nelem,1);

% Set Boundary Information
if meshtype ==2
    nsurf = 0;
    for i =1:nboun
        tempbrow = str2num(meshdata{9+i});
        j = tempbrow(1);
        k = tempbrow(2);
        i1 = tempbrow(3);
        i2 = tempbrow(4);
        np1 = tempbrow(5);
        bsido(i,1:3) = tempbrow(6:8);
        if j>nboun
            error(message('too many surfaces'));
        end
        if k==2 || k==3 || k==4 || k==5
            bcosurf(j) = k;
        else
            bcosurf(j) = 0;
        end
        if j >= nsurf
            nsurf = j;
        end
        if np1~=3 && np1~=6
            error(message('not a tetrahedral grid'));
        end
        bsido(i,5) = j;
        bsido(i,6)=k;
        bsido(i,4) = 0;
    end
    % detect automatically whether we have a linear mesh or quadratic
    % geometry
    if np1==3
        quadlin=1;
        disp('Mesh file has linear geometry')
    elseif np1==6
        quadlin=2;
        disp('Mesh file has quadratic geometry')
    end
    
else
    nsurf = 0;
    for i =1:nboun
        %        tempbrow = str2num(meshdata{9+i});
        tempbrow = str2num(meshdata{9+i});
        j = tempbrow(1);
        k = tempbrow(2);
        i1 = tempbrow(3);
        i2 = tempbrow(4);
        np1 = tempbrow(5);
        bsido(i,1:3) = tempbrow(6:8);
        if j>nboun
            error(message('too many surfaces'));
        end
        if k==2 || k==3 || k==4 || k==5
            bcosurf(j) = k;
        else
            bcosurf(j) = 0;
        end
        if j >= nsurf
            nsurf = j;
        end
        if np1~=3 && np1~=6
            error(message('not a tetrahedral grid'));
        end
        bsido(i,5) = j;
        bsido(i,6)=k;
        bsido(i,4) = 0;
    end
    % detect automatically whether we have a linear mesh or quadratic
    % geometry
    if np1==3
        quadlin=1;
        disp('Mesh file has linear geometry')
    elseif np1==6
        quadlin=2;
        disp('Mesh file has quadratic geometry')
    end
end

% Connectivities
nbounnew = nboun;
help1(1:nboun,1:6) = bsido(1:nboun,1:6);

for i =1:nelem
    tempcrow = str2num(meshdata{9+nboun+5+i});
    if quadlin ==1
        mat(i) = tempcrow(1);
        np2 = tempcrow(2);
        intma(i,1:4) = tempcrow(3:6);
    else
        mat(i) = tempcrow(1);
        np2 = tempcrow(2);
        intma(i,1:4) = tempcrow(3:6);
        cintma(i,1:6) = tempcrow(7:12);
    end
    if np2 ~= 4 && np2 ~=10
        error(message('not a tetrahedral mesh'));
    end
end

% add element flag to bsido
for i = 1:nboun
    ib = bsido(i,1:3);
    for j = 1:nelem
        jb = intma(j,1:4);
        match = 0;
        for k = 1:3
            for p = 1:4
                if ib(k) == jb(p)
                    match = match +1;
                end
            end
        end
        if match ==3
            if help1(i,4) ==0
                % display('element found');
                help1(i,4) =j;
            elseif help1(i,4) ~=0
                nbounnew = nbounnew +1;
                help1(nbounnew,4) =j;
                help1(nbounnew,1:3) = help1(i,1:3);
                help1(nbounnew,5) = help1(i,5);
                help1(nbounnew,6)=help1(i,6);
            end
        end
    end
    if help1(i,4)==0
        error(message(['No element found for this boundary face',num2str(i)]));
    end
end

%reset mat flag if necessary
if matflg ~=1
    for i = 1:nelem
        mat(i) = 1;
    end
end

if nbounnew > 2*nboun
    error(message('too many boundary faces'));
end

clear bsido
bsido = help1(1:nbounnew,1:6);


% Set up coordinates
ned = str2num(meshdata{9+nboun+5+nelem+5});
npoint = str2num(meshdata{9+nboun+5+nelem+5+ned+5});
coord = zeros(npoint,3);
for i=1:npoint
    coord(i,:) = str2num(meshdata{9+nboun+5+nelem+5+ned+5+i});
end
% apply shift and rescale
for i=1:npoint
     coord(i,1:3)=delta*(coord(i,1:3)-shift(1:3));
    %coord(i,1:3)=(coord(i,1:3)-shift(1:3));
    
end
nboun = nbounnew;
aa=intma(:,1);
bb=intma(:,2);
cc=intma(:,3);
dd=intma(:,4);
ee=[aa;bb;cc;dd];
nNodes=max(ee);






% Save relevant variables to mesh structure
mesh.nboun=nboun;
mesh.bsido=bsido;
mesh.Nelements=nelem;
mesh.intma=intma;
mesh.npoin=npoint;
mesh.nNodes=nNodes;
mesh.Coordinates=coord;
mesh.cintma=cintma;
mesh.mat=mat;
mesh.bcosurf=bcosurf;
mesh.nsurf=nsurf;
mesh.quadlin=quadlin;


