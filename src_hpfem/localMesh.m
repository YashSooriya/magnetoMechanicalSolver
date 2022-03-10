function [struct] = localMesh(mesh,probdata,physicType)
%=========================================================================
% Function to extract the local mesh data for the defined physic type from
% the global mesh data
%=========================================================================
% Inputs:
% mesh        : Data structure containing the generated global mesh data
% problemData : Data structure containing the problem data generated from
%               the problem file
% physicType  : The user-specified physic type for extracting the mesh
% plotsOn     : The plotting swith, used to determine if plotting occurs
% Outputs:
% struct      : Data structure containing the local mesh data

%-------------------------------------------------------------------------
% Extract the data from the inputs
%-------------------------------------------------------------------------
% Extract the global mesh data
nelem=mesh.Nelements;
nNodes=mesh.nNodes;
coord=mesh.Coordinates;
intma=mesh.intma;
cintma=mesh.cintma;
mat=mesh.mat;
matc=mesh.matc;
eltype=mesh.eltype;
nside=mesh.edge.nside;
glob=mesh.edge.glob;
nface=mesh.face.nface;
globfa=mesh.face.globfa;
bcedge=mesh.edge.bcedge;
bcvertex=mesh.bcvertex;
cond=mesh.face.cond;
subFlag=mesh.subFlag;
[dum,cint_size]=size(cintma);
quadlin=mesh.quadlin;



% Determine the number of local mesh elements
nElemNew    = sum(subFlag==physicType);


%-------------------------------------------------------------------------
% Compute a set of mapping variables for the local to global mesh
%-------------------------------------------------------------------------
% Initialise elemental variables
nintma     = zeros(nElemNew,4);           % Triangular Elements
subdomNew     = zeros(nElemNew,1);

mapG2L_e    = zeros(1,nelem);       % Global Mappings
mapL2G_e    = zeros(1,nElemNew);

% Initialise loop counters
elem        = 0;


% Initialise nodal variables
help        = zeros(length(coord),1);
ncoord      = [];
mapG2L_n    = zeros(1,length(coord));
mapL2G_n    = zeros(1,1);
node        = 0;

% Initialise edge variables
helpedge    = zeros(nside,1);
helpface    = zeros(nface,1);

nglob       = zeros(nElemNew,6);
nglobfa     = zeros(nElemNew,4);
mapG2L_ed   = zeros(1,nside);
mapL2G_ed   = zeros(1,1);
edge        = 0;
face        = 0;

for i=1:nelem
    % Determine if the element lies in the mechanical subdomain
    if subFlag(i)==physicType
        %-----------------------------------------------------------------
        % Map the global to local elements
        %-----------------------------------------------------------------
        % Determine the submesh element number 
        elem = elem+1;
        eltype_mech(elem)=eltype(i);
        matc_mech(elem)=matc(i);
        mat_mech(elem)=mat(i);
        
        % Store the global to local element number
        mapG2L_e(i) = elem;
        
        % Store the local to global element number
        mapL2G_e(elem) = i;
        
        %-----------------------------------------------------------------
        % Loop over edges and nodes
        %-----------------------------------------------------------------
        % Initialise the edge variable
        e = zeros(4,1);
        
    
            temp  = intma(i,:);
            temp2 = glob(i,:);
            temp3 = globfa(i,:);
            if quadlin>1
            tempc = cintma(i,:);
            end
            npel  = 4;
            edpel = 6;
            fapel = 4;
            c     = zeros(1,npel);
            
      % Boundary information
                condNew(elem,:)=cond(i,:);

        % Loop through the nodes
        for j=1:npel
            %-------------------------------------------------------------
            % Map the global to local nodes
            %-------------------------------------------------------------
            if help(temp(j))==0
                % Determine the local node number
                node = node+1;
                
                % Store the local node information
                c(j) = node;
                
                % Store the nodal coordinate information
                ncoord(node,:) = coord(temp(j),:);
               bcvertexn(node)=bcvertex(temp(j));
                
                % Update the node connectivity identifier
                help(temp(j)) = node;
                
                % Store the global to local element number
                mapG2L_n(temp(j)) = node;
                
                % Store the local to global element number
                mapL2G_n(node) = temp(j);
                
                
                
            else
                % Store the local node information
                c(j) = help(temp(j));
            end
        end

        
        for j=1:edpel
            
            %-------------------------------------------------------------
            % Map the global to local edges
            %-------------------------------------------------------------
            if helpedge(temp2(j))==0
                % Determine the local edge number
                edge = edge + 1;
                
                % Store the local edge information
                e(j) = edge;
                
                % Bouncary condition information
                bcedgeNew(edge)=bcedge(temp2(j));
                % Update the node connectivity identifier
                helpedge(temp2(j)) = edge;
                
                % Store the global to local element number
                mapG2L_ed(temp2(j)) = edge;
                
                % Store the local to global element number
                mapL2G_ed(edge) = temp2(j);
                
            else
                % Store the local edge information
                e(j) = helpedge(temp2(j));
            end
            
        end
        
        for j=1:fapel
                if helpface(temp3(j))==0
                % Determine the local edge number
                face = face + 1;
                
                % Store the local edge information
                f(j) = face;
                
            
                
                % Update the node connectivity identifier
                helpface(temp3(j)) = face;
                
                % Store the global to local element number
                mapG2L_fa(temp3(j)) = face;
                
                % Store the local to global element number
                mapL2G_fa(face) = temp3(j);
                
                else
                % Store the local edge information
                f(j) = helpface(temp3(j));
                end
            
        end

        %-----------------------------------------------------------------
        % Store the elemental data 
        %-----------------------------------------------------------------
   
            % Store the sub domain flag
            subdomNew(elem) = subFlag(i);
            
            % Store data in the local connectivity array
            nintma(elem,:) = c;
       
                
        % store data in the local edge array
        nglob(elem,:) = e;
        nglobfa(elem,:)=f;
    end
end

helpc=zeros(length(coord),1);
  if quadlin>1 
       nodeQuad=node;
  k=0;
for i=1:nelem
    if subFlag(i)==physicType
       k=k+1;
         tempc = cintma(i,:);
           temp2 = glob(i,:);
            temp3 = globfa(i,:);
 
 for j=1:cint_size
     if helpc(tempc(j))==0
     nodeQuad=nodeQuad+1;
        cintma_mech(k,j)=nodeQuad;
        ncoord(nodeQuad,:)=coord(tempc(j),:);
        helpc(tempc(j))=nodeQuad;
        mapG2L_n(tempc(j))=nodeQuad;
        mapL2G_n(nodeQuad)=tempc(j);
     else
         cintma_mech(k,j)=helpc(tempc(j));
     end
 end
   end
end
  end


aa=nintma(:,1);
bb=nintma(:,2);
cc=nintma(:,3);
dd=nintma(:,4);
ee=[aa;bb;cc;dd];
%nNodes=max(mapL2G_n(ee));
nNodes=node;




%-------------------------------------------------------------------------
% Store the mechanical mesh data and plot the mesh
%-------------------------------------------------------------------------
struct.mapG2L_e=mapG2L_e;
struct.mapL2G_e=mapL2G_e;
struct.Coordinates=ncoord;
struct.mapG2L_n=mapG2L_n;
struct.mapL2G_n=mapL2G_n;
struct.bcedge=bcedgeNew;
struct.mapG2L_ed=mapG2L_ed;
struct.mapL2G_ed=mapL2G_ed;
struct.cond=condNew;
struct.mapG2L_fa=mapG2L_fa;
struct.mapL2G_fa=mapL2G_fa;
struct.subdom=subdomNew;
struct.intma=nintma;
struct.nface=face;
struct.nside=edge;
struct.glob=nglob;
struct.globfa=nglobfa;
struct.Nelements=elem;
struct.npoin=node;
struct.bcvertex=bcvertexn;
struct.nNodes=nNodes;
struct.eltype=eltype_mech;
struct.mat=mat_mech;
struct.matc=matc_mech;
struct.face.nface=face;
if quadlin>1
struct.cintma=cintma_mech;
end
end