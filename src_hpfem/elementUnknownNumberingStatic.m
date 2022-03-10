function [Unknown]=elementUnknownNumberingStatic(mesh,ProblemData)

%-------------------------------------------------------------------------
% Extract the data from de data structures
%-------------------------------------------------------------------------
nelem=mesh.Nelements;
nelemM=mesh.mech.Nelements;
order=ProblemData.order;
orderH1=ProblemData.orderH1;
esizet=ProblemData.esizet;
esizeH1=ProblemData.esizeH1;
subFlag=mesh.subFlag;


% Estimate the size of sparse matrix vectors
nSparse=floor((((nelem-nelemM)*esizet^2+nelemM*(esizet+2.25*esizeH1)^2)));
%nSparse=nelemM*(3*esizeH1)^2+nelem*(esizet)^2;

%-------------------------------------------------------------------------
% Determine the global unknown numbering
%-------------------------------------------------------------------------
if ProblemData.probstatic==1
    [Unknown]=nounk_static(mesh,ProblemData);
else
    [Unknown]=nounk_static3(mesh,ProblemData);
end

%-------------------------------------------------------------------------
% Determine the elemental unknown numbering
%-------------------------------------------------------------------------
% Initialise the elemental unknown numbering matrices
lunkvEM=zeros(esizet,nelem);
lunkvM=zeros(3*esizeH1,nelem);


% Loop over elements
for i=1:nelem
    [lunkvEM]=elementalUnk(Unknown,order,orderH1,mesh,i,subFlag(i),lunkvEM,1);
    [lunkvM]=elementalUnk(Unknown,order,orderH1,mesh,i,subFlag(i),lunkvM,2);
end

%-------------------------------------------------------------------------
% Adjust the unknowns so that Dirichlet numbers are corrected
%-------------------------------------------------------------------------

% Adjust the electromagnetic system
lunkvEMf=lunkvEM;
lunkvEMf(lunkvEMf<0)=abs(lunkvEMf(lunkvEMf<0))+Unknown.EM.nunkt;

% Adjust the mechanical system
lunkvMf=lunkvM;
lunkvMf(lunkvMf<0)=abs(lunkvMf(lunkvMf<0))+Unknown.Mech.nunkt;


%-------------------------------------------------------------------------
% Store data in the Unknown structure
%-------------------------------------------------------------------------

% Store the electromagnetic unknowns
Unknown.EM.unknowns=lunkvEMf;

% Store the mechanical unknowns
Unknown.Mech.unknowns=lunkvMf;


% Correct the mechanical numbering to store in system
lunkvM(lunkvM>0)=lunkvM(lunkvM>0)+Unknown.EM.nunkt;
lunkvM(lunkvM<0)=lunkvM(lunkvM<0)-Unknown.EM.npec;

% Store the system unknowns
Unknown.system.unknowns=[lunkvEM;lunkvM];

% Correct the mechanical system Dirichlet values
lunkvMs=lunkvM;
lunkvMs(lunkvMs>0)=lunkvMs(lunkvMs>0)+Unknown.EM.npec;
lunkvMs(lunkvMs<0)=abs(lunkvMs(lunkvMs<0))+Unknown.EM.nunkt+Unknown.Mech.nunkt;
Unknown.Mech.unknownsD=lunkvMs;

% Store size of sparse matrix vectors
Unknown.system.nSparse=nSparse;

%=========================================================================
function out=elementalUnk(unknown,order,orderH1,mesh,i,subFlag,out,probFlag)
%=========================================================================
% Determine the element unknown numbers
%-------------------------------------------------------------------------
% Electromagntic DOF numbering
%-------------------------------------------------------------------------
if probFlag==1
    
    % Extract relevant data from structures
    unkz=unknown.EM.unkz;
    unksid=unknown.EM.unksid;
    unkint=unknown.EM.unkint;
    unkfatp1=unknown.EM.unkfatp1;
    unkfatp2=unknown.EM.unkfatp2;
    unkfatp3=unknown.EM.unkfatp3;
    
    glob=mesh.edge.glob;
    globfa=mesh.face.globfa;
    
    % Initialise counting variable
    nbas=0;
    
    % Low order edge basis
    for m=1:6
        nbas=nbas+1;
        out(nbas,i)=unkz(glob(i,m));
    end
    
    % Higher order edge basis
    if order>=1
        for pp=1:order
            for m=1:6
                nbas=nbas+1;
                out(nbas,i)=unksid(glob(i,m),pp);
            end
        end
    end
    
    % Face basis
    if order>=2
        
        for m=1:4
            
            for k=1:(order*order-order)/2
                nbas=nbas+1;
                out(nbas,i)=unkfatp1(globfa(i,m),k);
            end
        end
        
        for m=1:4
            
            for k=1:(order*order-order)/2
                nbas=nbas+1;
                out(nbas,i)=unkfatp2(globfa(i,m),k);
            end
        end
        for m=1:4
            for jj=0:order-2
                nbas=nbas+1;
                out(nbas,i)=unkfatp3(globfa(i,m),jj+1);
            end
        end
        
    end
    
    % Interior functions
    if order>= 3
        k=0;
        for ii=0:order-3
            for jj=0:order-3
                for kk=0:order-3
                    if ii+jj+kk <= order-3
                        nbas=nbas+1;
                        k=k+1;
                        out(nbas,i)=unkint(i,k);
                    end
                end
            end
        end
        
        % Non-gradients
        nintbas=k;
        for ii=nintbas+1:(order-2)*(order-1)*(order+1)/2
            k=k+1;
            nbas=nbas+1;
            out(nbas,i)=unkint(i,k);
            
        end
    end
end
%-------------------------------------------------------------------------
% Inlude mechanical (H1) degrees of freedom
%-------------------------------------------------------------------------
if probFlag==2
    %if subFlag==2
    
    % Extract relevant data from structures
    unkvertexMech=unknown.Mech.unkvertexMech;
    unkedgesMech=unknown.Mech.unkedgesMech;
    unkfacesMech=unknown.Mech.unkfacesMech;
    unkinteriorsMech=unknown.Mech.unkinteriorsMech;
    unkvertexSystem=unknown.system.unkvertexSystem;
    unkedgesSystem=unknown.system.unkedgesSystem;
    unkfacesSystem=unknown.system.unkfacesSystem;
    unkinteriorsSystem=unknown.system.unkinteriorsSystem;
    intmaMech=mesh.mech.intma;
    intma=mesh.intma;
    glob=mesh.mech.glob;
    globTotal=mesh.edge.glob;
    globfaTotal=mesh.face.globfa;
    globfa=mesh.mech.globfa;
    mapG2L_e=mesh.mech.mapG2L_e;
    mapL2G_n=mesh.mech.mapL2G_n;
    
    % Mapping global to local (mechanical) elements
    iii=i;
    i=mapG2L_e(i);
    
    % Initialise counting variable
    nbas=0;
    
    % Vertex basis
    for m=1:4
        for k=1:3
            nbas=nbas+1;
            
            %out(nbas)=unkvertexSystem((intma(iii,m)),k)-unknown.EM.nunkt;
            out(nbas,iii)=unkvertexMech(intma(iii,m),k);
        end
    end
    
    % Edge basis
    
    if orderH1>=2
        
        for j=0:orderH1-2
            for m=1:6
                for k=1:3
                    nbas=nbas+1;
                    %out(nbas,iii)=unkedgesSystem(globTotal(iii,m),j+1,k)-unknown.EM.nunkt;
                    out(nbas,iii)=unkedgesMech(globTotal(iii,m),j+1,k);
                end
            end
        end
    end
    
    % Face basis
    
    if orderH1>=3
        for m=1:4
            for j=1:((orderH1-1)^2-(orderH1-1))/2
                for k=1:3
                    nbas=nbas+1;
                    %out(nbas,iii)=unkfacesSystem(globfaTotal(iii,m),j,k)-unknown.EM.nunkt;
                    out(nbas,iii)=unkfacesMech(globfaTotal(iii,m),j,k);
                end
            end
        end
    end
    
    %Interior basis
    
    if orderH1>=4
        k=0;
        for ii=0:orderH1-4
            for jj=0:orderH1-4
                for kk=0:orderH1-4
                    if ii+jj+kk <= orderH1-4
                        k=k+1;
                        for l=1:3
                            nbas=nbas+1;
                            % out(nbas,iii)=unkinteriorsSystem(iii,k,l)-unknown.EM.nunkt;
                            out(nbas,iii)=unkinteriorsMech(iii,k,l);
                        end
                    end
                end
            end
        end
    end
    % end
end