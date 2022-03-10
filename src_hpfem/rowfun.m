function out=rowfun(unknown,order,orderH1,mesh,i,esize,esizeH1,subFlag)

unkz=unknown.EM.unkz;
unksid=unknown.EM.unksid;
unkint=unknown.EM.unkint;
unkfatp1=unknown.EM.unkfatp1;
unkfatp2=unknown.EM.unkfatp2;
unkfatp3=unknown.EM.unkfatp3;
nunkt=unknown.EM.nunkt;

glob=mesh.edge.glob;
globfa=mesh.face.globfa;
                    
                    
nbas=0;
out = zeros(esize+3*esizeH1,1);
% edge basis functions
for m=1:6
    nbas=nbas+1;
    out(nbas)=unkz(glob(i,m));
end

if order>=1
    for pp=1:order
        for m=1:6
            nbas=nbas+1;
            out(nbas)=unksid(glob(i,m),pp);
        end
    end
end

% face basis functions
if order>=2
    
    for m=1:4

        for k=1:(order*order-order)/2
              nbas=nbas+1;
              out(nbas)=unkfatp1(globfa(i,m),k);
        end
    end
    
    for m=1:4

        for k=1:(order*order-order)/2
              nbas=nbas+1;
              out(nbas)=unkfatp2(globfa(i,m),k);
        end
    end
    for m=1:4
        for jj=0:order-2
            nbas=nbas+1;
            out(nbas)=unkfatp3(globfa(i,m),jj+1);
        end
    end

end
if order>= 3
% interior functions
    k=0;
    for ii=0:order-3
        for jj=0:order-3
            for kk=0:order-3
                if ii+jj+kk <= order-3
                    nbas=nbas+1;
                    k=k+1;
                    out(nbas)=unkint(i,k);
                end
            end
        end
    end
    
    %     non-gradients
    nintbas=k;
    for ii=nintbas+1:(order-2)*(order-1)*(order+1)/2
        k=k+1;
        nbas=nbas+1;
        out(nbas)=unkint(i,k);
        
    end
end

%-----------------------------------------------------------------------------------------------
% Inlude mechanical (H1) degrees of freedom
%-----------------------------------------------------------------------------------------------
if subFlag==2
    
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
globfa=mesh.mech.globfa;
mapG2L_e=mesh.mech.mapG2L_e;
mapL2G_n=mesh.mech.mapL2G_n;

  i=mapG2L_e(i);

% Vertex
for m=1:4
    for k=1:3
    nbas=nbas+1;
    
    %out(nbas)=unkvertexSystem(mapL2G_n(intmaMech(i,m)),k);
    out(nbas)=unkvertexSystem(intmaMech(i,m),k);
    end
end

% Edge

if orderH1>=2

    for j=0:orderH1-2
        for m=1:6
            for k=1:3
    nbas=nbas+1;
    out(nbas)=unkedgesSystem(glob(i,m),j+1,k);
            end
        end
    end
end

% Face

if orderH1>=3
    for m=1:4
        for j=1:((orderH1-1)^2-(orderH1-1))/2
            for k=1:3
        nbas=nbas+1;
        out(nbas)=unkfacesSystem(globfa(i,m),j,k);
            end
        end
    end
end

%Interiors

if orderH1>=4
    k=0;
            for ii=0:orderH1-4
                for jj=0:orderH1-4
                    for kk=0:orderH1-4
                        if ii+jj+kk <= orderH1-4
                            k=k+1;
                            for l=1:3
                            nbas=nbas+1;
                            out(nbas)=unkinteriorsSystem(i,k,l);
                            end
                        end
                    end
                end
            end
end
end