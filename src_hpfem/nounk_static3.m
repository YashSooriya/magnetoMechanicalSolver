function [unknown]=nounk_static3(mesh,ProblemData)

% Extract relevant data from mesh structure
nelem=mesh.Nelements;
nside=mesh.edge.nside;
nface=mesh.face.nface;
cond=mesh.face.cond;
glob=mesh.edge.glob;
globfa=mesh.face.globfa;
bcedge=mesh.edge.bcedge;
matc=mesh.matc;
matc=0*matc;


% Extract data from FEspacesInfo structure
order=ProblemData.order;

% unksid = unknown numbers on element edges
% unkint = unknown numbers on element interiors
% unkfac = unknown numbers on element faces

% set up help array
help1 = zeros(1,nface);
for i=1:nelem
    for j=1:4
        if  cond(i,j)==2 || cond(i,j)==10 || cond(i,j)==12
            help1(globfa(i,j))=cond(i,j);
        end
    end
end


% flag edges and faces in Region 1 so as to skip gradients
help2 = zeros(1,nside);
help3 = zeros(1,nface);

% Flag only those edges in the interior of free space
for i=1:nelem
    if matc(i)==0
        for j=1:6
            help2(glob(i,j))=1;
        end
    end
end
for i=1:nelem
    if matc(i)==1
        for j=1:6
            help2(glob(i,j))=0;
        end
    end
end

% Flag only those faces in the interior of free space
for i=1:nelem
    if matc(i)==0
        for j=1:4
            help3(globfa(i,j))=1;
        end
    end
end
for i=1:nelem
    if matc(i)==1
        for j=1:4
            help3(globfa(i,j))=0;
        end
    end
end




unkz = zeros(nside,1);
unksid = zeros(nside,order);

npec=0;
nunk=0;
nintbas=0;


% treat 1,2,4 Dirichlet (2 posibly non-zero)
% treat 3 Neumann
for i = 1:nside
    if  bcedge(i)~=2 && bcedge(i)~=10 && bcedge(i)~=12 && bcedge(i)~=612  && bcedge(i)~=128  && bcedge(i)~=610
        nunk = nunk +1;
        unkz(i) = nunk;
    elseif bcedge(i) ==2 || bcedge(i)==10 || bcedge(i)==12 || bcedge(i)==612 || bcedge(i)==128 || bcedge(i)==610 
        npec = npec+1;
        unkz(i) = -1*npec;
    end
end

% higher order edges
if order >=1
    for i=1:nside
        %        if help2(i)==0
        if bcedge(i)~=2 && bcedge(i)~=10 && bcedge(i)~=12 && bcedge(i)~=612  && bcedge(i)~=128  && bcedge(i)~=610  && help2(i)==0
            for j = 1:order
                nunk = nunk+1;
                unksid(i,j) = nunk;
            end 
        elseif bcedge(i) ==2 || bcedge(i)==10 || bcedge(i)==12 || bcedge(i)==612 || bcedge(i)==128 || bcedge(i)==610 
            for j=1:order
                npec = npec+1;
                unksid(i,j) = -1*npec;
            end
        end
    end
    %    end
end
display(['The number of edge based unknowns is ',num2str(nunk)])

% Face functions
unkfatp1=[];
unkfatp2=[];
unkfatp3=[];
if order >=2
    % zero unkfatp1
    unkfatp1 = zeros(nface,(order*order-order)/2);
    for i =1:nface
        %       if help3(i)==0
        if help1(i) ==0 && help3(i)==0
            k = 0;
            for j1 = 0:order-2
                for j2=0:order-2
                    if (j1+j2)<=(order-2)
                        k=k+1;
                        nunk = nunk+1;
                        unkfatp1(i,k)=nunk;
                    end
                end
            end
        elseif help1(i)==2 || help1(i)==10 || help1(i)==12
            k=0;
            for j1=0:order-2
                for j2=0:order-2
                    if (j1+j2)<=(order-2)
                        k = k+1;
                        npec=npec+1;
                        unkfatp1(i,k)=-1*npec;
                    end
                end
            end
            
        end
        %      end
    end
    
    % zero unkfatp2
    unkfatp2 = zeros(nface,(order*order-order)/2);
    
    % zero unkfatp3
    unkfatp3 = zeros(nface,order-1);
    for i=1:nface
        if help1(i)==0
            k=0;
            for j1=0:order-2
                for j2=0:order-2
                    if (j1+j2)<=(order-2)
                        k=k+1;
                        nunk = nunk+1;
                        unkfatp2(i,k)=nunk;
                    end
                end
            end
        elseif help1(i)==2 || help1(i)==10 || help1(i)==12
            k=0;
            for j1=0:order-2
                for j2=0:order-2
                    if (j1+j2)<=(order-2)
                        k=k+1;
                        npec=npec+1;
                        unkfatp2(i,k)=-1*npec;
                    end
                end
            end
        end
        %end
        
        %for i = 1:nface
        if help1(i)==0
            for j=1:order-1
                nunk=nunk+1;
                unkfatp3(i,j)= nunk;
            end
        elseif help1(i)==2 || help1(i)==10 || help1(i)==12
            for j = 1:order-1
                npec = npec+1;
                unkfatp3(i,j)=-1*npec;
            end
        end
    end
end
display(['The total number of edge and face unknowns is',num2str(nunk)]);
display(['The total number of PEC edge and face knowns is',num2str(npec)]);

% save number of H curl unknowns
nef = nunk;
nunkt = nunk;

% Interior functions
unkint = [];
if order >=3
    unkint = zeros(nelem,(order-1)*(order-2)*(order+1)/2);
    for i = 1:nelem
        if matc(i)==1
            % gradients
            k=0;
            for ii=0:order-3
                for jj=0:order-3
                    for kk=0:order-3
                        if ii+jj+kk <= order-3
                            k=k+1;
                            nunkt = nunkt+1;
                            unkint(i,k) = nunkt;
                        end
                    end
                end
            end
        end
    end
    % if there are no conductors the above loop may give k=0
    % so as to fix the numbering do the following
    k=0;
    for ii=0:order-3
        for jj=0:order-3
            for kk=0:order-3
                if ii+jj+kk <= order-3
                    k=k+1;
                end
            end
        end
    end
    % non -gradients (do not distinguish between)
    nintbas=k;
    for i = 1:nelem
        k=nintbas;
        for ii=nintbas+1:(order-2)*(order-1)*(order+1)/2
            nunkt = nunkt+1;
            k=k+1;
            unkint(i,k) = nunkt;
        end
    end
    
    
    if k~=(order-2)*(order-1)*(order+1)/2
        disp('wrong number of interior functions created')
    end
end

% Save relevant data to unknown structure
unknown.EM.nunk=nunk;
unknown.EM.nunkt=nunkt;
unknown.EM.unkz=unkz;
unknown.EM.unksid=unksid;
unknown.EM.unkfatp1=unkfatp1;
unknown.EM.unkfatp2=unkfatp2;
unknown.EM.unkfatp3=unkfatp3;
unknown.EM.unkint=unkint;
unknown.EM.nef=nef;
unknown.EM.npec=npec;
unknown.EM.nintbas=nintbas;

display(['The total number of edge, face and interior unknowns is',num2str(nunkt)]);

%=========================== Include H1 degrees of freedom=============

% Extract data from mechanical mesh structure
nelem=mesh.Nelements;
nside=mesh.edge.nside;
nface=mesh.face.nface;
cond=mesh.face.cond;
glob=mesh.edge.glob;
globfa=mesh.face.globfa;
bcedge=mesh.edge.bcedge;
intma=mesh.intma;
bcvertex=mesh.bcvertex;
nNodes=mesh.nNodes;
FlagMechVolume=zeros(nelem,1);


% Define order for H1 space
orderH1=ProblemData.orderH1;

% Vertex unknowns

nunkMech=0;
npecMech=0;

subFlag=mesh.subFlag;


unkvertexMech=zeros(nNodes,3);
unkvertexSystem=zeros(nNodes,3);
for i=1:nNodes
    [row,~]=find(intma==i);
    if any(subFlag(row)==2)
        FlagMechVolume(row)=1;
        % ii=mapG2L_n(i);
        %if ii~=0
        if bcvertex(i)~=6 && bcvertex(i)~=7 && bcvertex(i)~=8
            for k=1:3
                nunkMech=nunkMech+1;
                nunkt=nunkt+1;
                unkvertexMech(i,k)=nunkMech;
                unkvertexSystem(i,k)=nunkt;
            end
        elseif bcvertex(i)==6
            for k=1:3
                npecMech = npecMech+1;
                npec=npec+1;
                unkvertexMech(i,k) = -1*npecMech;
                unkvertexSystem(i,k)=-1*npec;
            end
            
        elseif bcvertex(i)==7
            nunkMech=nunkMech+1;
            nunkt=nunkt+1;
            unkvertexMech(i,2)=nunkMech;
            unkvertexSystem(i,2)=nunkt;
            nunkMech=nunkMech+1;
            nunkt=nunkt+1;
            unkvertexMech(i,3)=nunkMech;
            unkvertexSystem(i,3)=nunkt;
            npecMech=npecMech+1;
            npec=npec+1;
            unkvertexMech(i,1)=-1*npecMech;
            unkvertexSystem(i,1)=-1*npec;
        elseif bcvertex(i)==8
            nunkMech=nunkMech+1;
            nunkt=nunkt+1;
            unkvertexMech(i,1)=nunkMech;
            unkvertexSystem(i,1)=nunkt;
            nunkMech=nunkMech+1;
            nunkt=nunkt+1;
            unkvertexMech(i,3)=nunkMech;
            unkvertexSystem(i,3)=nunkt;
            npecMech=npecMech+1;
            npec=npec+1;
            unkvertexMech(i,2)=-1*npecMech;
            unkvertexSystem(i,2)=-1*npec;
        end
        
    end
end

nunkVertex=nunkt;
% Edge unknowns
help2=zeros(1,nside);




unkedgesMech=[];
unkedgesSystem=[];
if orderH1>=2
    unkedgesMech=zeros(nside,orderH1-1,3);
    unkedgesSystem=zeros(nside,orderH1-1,3);
    
    for i=1:nside
        [row,~]=find(glob==i);
        if any(subFlag(row)==2)
            FlagMechVolume(row)=1;
            for j=0:orderH1-2
                
                
                if bcedge(i)~=6 && bcedge(i)~=7 && bcedge(i)~=8
                    for k=1:3
                        nunkMech=nunkMech+1;
                        nunkt=nunkt+1;
                        unkedgesMech(i,j+1,k)=nunkMech;
                        unkedgesSystem(i,j+1,k)=nunkt;
                    end
                elseif bcedge(i)==6
                    for k=1:3
                        npecMech=npecMech+1;
                        npec=npec+1;
                        unkedgesMech(i,j+1,k)=-1*npecMech;
                        unkedgesSystem(i,j+1,k)=-1*npec;
                    end
                elseif bcedge(i)==7
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkedgesMech(i,j+1,2)=nunkMech;
                    unkedgesSystem(i,j+1,2)=nunkt;
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkedgesMech(i,j+1,3)=nunkMech;
                    unkedgesSystem(i,j+1,3)=nunkt;
                    npecMech=npecMech+1;
                    npec=npec+1;
                    unkedgesMech(i,j+1,1)=-1*npecMech;
                    unkedgesSystem(i,j+1,1)=-1*npec;
                elseif bcedge(i)==8
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkedgesMech(i,j+1,1)=nunkMech;
                    unkedgesSystem(i,j+1,1)=nunkt;
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkedgesMech(i,j+1,3)=nunkMech;
                    unkedgesSystem(i,j+1,3)=nunkt;
                    npecMech=npecMech+1;
                    npec=npec+1;
                    unkedgesMech(i,j+1,2)=-1*npecMech;
                    unkedgesSystem(i,j+1,2)=-1*npec;
                    
                end
                
            end
        end
    end
end
nunkEdges=nunkt;

% Face unknowns

% Create help array so that it is zero only if cond(i,j)=0 or
% cond(i,j)=3
helpF=zeros(1,nface);


for i=1:nelem
    for j=1:4
        if  cond(i,j)==6 || cond(i,j)==7 || cond(i,j)==8
            helpF(globfa(i,j))=cond(i,j);
        end
    end
end



unkfacesMech=[];
unkfacesSystem=[];
if orderH1>=3
    unkfacesMech = zeros(nface,((orderH1-1)^2-(orderH1-1))/2,3);
    unkfacesSystem = zeros(nface,((orderH1-1)^2-(orderH1-1))/2,3);
    for i =1:nface
        [row,~]=find(globfa==i);
        if any(subFlag(row)==2)
            FlagMechVolume(row)=1;
            for j=1:((orderH1-1)^2-(orderH1-1))/2
                if  helpF(i)==0
                    
                    for k=1:3
                        nunkMech=nunkMech+1;
                        nunkt=nunkt+1;
                        unkfacesMech(i,j,k)=nunkMech;
                        unkfacesSystem(i,j,k)=nunkt;
                        
                    end
                    
                elseif helpF(i)==6
                    for k=1:3
                        npecMech=npecMech+1;
                        npec=npec+1;
                        unkfacesMech(i,j,k)=-1*npecMech;
                        unkfacesSystem(i,j,k)=-1*npec;
                    end
                elseif helpF(i)==7
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkfacesMech(i,j,2)=nunkMech;
                    unkfacesSystem(i,j,2)=nunkt;
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkfacesMech(i,j,3)=nunkMech;
                    unkfacesSystem(i,j,3)=nunkt;
                    npecMech=npecMech+1;
                    npec=npec+1;
                    unkfacesMech(i,j,1)=-1*npecMech;
                    unkfacesSystem(i,j,1)=-1*npec;
                elseif helpF(i)==8
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkfacesMech(i,j,1)=nunkMech;
                    unkfacesSystem(i,j,1)=nunkt;
                    nunkMech=nunkMech+1;
                    nunkt=nunkt+1;
                    unkfacesMech(i,j,3)=nunkMech;
                    unkfacesSystem(i,j,3)=nunkt;
                    npecMech=npecMech+1;
                    npec=npec+1;
                    unkfacesMech(i,j,2)=-1*npecMech;
                    unkfacesSystem(i,j,2)=-1*npec;
                    
                end
            end
        end
    end
end
nunkFaces=nunkt;
% Interior unknowns

unkinteriorsMech=[];
unkinteriorsSystem=[];
if orderH1>=4
    unkinteriorsMech=zeros(nelem, (orderH1-3)*(orderH1-2)*(orderH1-1)/6,3);
    unkinteriorsSystem=zeros(nelem, (orderH1-3)*(orderH1-2)*(orderH1-1)/6,3);
    for i = 1:nelem
        if subFlag(i)==2
            FlagMechVolume(i)=1;
            k=0;
            for ii=0:orderH1-4
                for jj=0:orderH1-4
                    for kk=0:orderH1-4
                        if ii+jj+kk <= orderH1-4
                            k=k+1;
                            for l=1:3
                                nunkMech=nunkMech+1;
                                nunkt=nunkt+1;
                                unkinteriorsMech(i,k,l)=nunkMech;
                                unkinteriorsSystem(i,k,l)=nunkt;
                                
                            end
                        end
                    end
                end
            end
        end
    end
end
nunkInteriors=nunkt;

unknown.Mech.nunkt=nunkMech;
unknown.Mech.unkvertexMech=unkvertexMech;
unknown.Mech.unkedgesMech=unkedgesMech;
unknown.Mech.unkfacesMech=unkfacesMech;
unknown.Mech.unkinteriorsMech=unkinteriorsMech;
unknown.Mech.npec=npecMech;
unknown.Mech.nunkVertex=nunkVertex;
unknown.Mech.nunkEdges=nunkEdges;
unknown.Mech.nunkFaces=nunkFaces;
unknown.Mech.nunkInteriors=nunkInteriors;

unknown.system.nunkt=nunkt;
unknown.system.unkvertexSystem=unkvertexSystem;
unknown.system.unkedgesSystem=unkedgesSystem;
unknown.system.unkfacesSystem=unkfacesSystem;
unknown.system.unkinteriorsSystem=unkinteriorsSystem;
unknown.FlagMechVolume=FlagMechVolume;


display(['The total number of EM and mechanical unknowns is ',num2str(nunkt)]);

