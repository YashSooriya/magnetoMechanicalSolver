function [stiffVV,stiffEE,stiffFF,stiffII,nVV,nEE,nFF,nII]=extract_mech2(unknown,K,ProblemData)

%=========================================================================
% Extract data from structures
%=========================================================================
orderH1=ProblemData.orderH1;
nunkVertex=unknown.Mech.nunkVertex;
nunkEdges=unknown.Mech.nunkEdges;
nunkFaces=unknown.Mech.nunkFaces;
nunkInteriors=unknown.Mech.nunkInteriors;
nunkt=unknown.EM.nunkt;
unkEdgesMech=unknown.system.unkedgesSystem;
unkFacesMech=unknown.system.unkfacesSystem;
unkInteriorsMech=unknown.system.unkinteriorsSystem;

% Initialise matrices
stiffEE=[];
stiffFF=[];
stiffII=[];
nEE=0;
nFF=0;
nII=0;


%=========================================================================
% Mechanical part
%=========================================================================

nstart=1;
nend=nunkVertex-nunkt;

stiffVV=K(nstart:nend,nstart:nend);
nVV=nend;



if orderH1>1
    nstart=nend+1;
    nend=nunkEdges-nunkt;
    nEE=nend-nVV;
    
    stiffEE=K(nstart:nend,nstart:nend);
    [stiffEE]=extractsparse(stiffEE,unkEdgesMech,nstart-1+nunkt);
end

if orderH1>2
    nstart=nend+1;
    nend=nunkFaces-nunkt;
    nFF=nend-nVV-nEE;
    stiffFF=K(nstart:nend,nstart:nend);
end
[stiffFF]=extractsparse(stiffFF,unkFacesMech,nstart-1+nunkt);


if orderH1>3
    nstart=nend+1;
    nend=nunkInteriors-nunkt;
    stiffII=K(nstart:nend,nstart:nend);
    nII=nend-nVV-nFF-nEE;
    [stiffII]=extractsparse(stiffII,unkInteriorsMech,nstart-1+nunkt);
end








