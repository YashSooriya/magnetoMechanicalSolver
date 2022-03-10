function [unknown,sol]=streat(mesh,unknown,known,known_mech,sol,FEspacesInfo)

% Extract relevant data from unknown structure
npec=unknown.EM.npec;
nunkt=unknown.EM.nunkt;
unksid=unknown.EM.unksid;
unkfatp1=unknown.EM.unkfatp1;
unkfatp2=unknown.EM.unkfatp2;
unkfatp3=unknown.EM.unkfatp3;
unkz=unknown.EM.unkz;
unkint=unknown.EM.unkint;
nunktMech=unknown.Mech.nunkt;

% Extract relevant data from mesh structure
nside=mesh.edge.nside;
nface=mesh.face.nface;
nelem=mesh.Nelements;

% Extract relevant data from FEspacesInfo structure
order=FEspacesInfo.order;
orderH1=FEspacesInfo.orderH1;

for i=1:nside
    if unkz(i)<0
        unkz(i)=abs(unkz(i))+nunkt+nunktMech;
    end
end

for i=1:nside
    for j=1:order
        if unksid(i,j)<0
            unksid(i,j)=abs(unksid(i,j))+nunkt+nunktMech;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if unkfatp1(i,j)<0
            unkfatp1(i,j)=abs(unkfatp1(i,j))+nunkt+nunktMech;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if(unkfatp2(i,j)<0)
            unkfatp2(i,j)=abs(unkfatp2(i,j))+nunkt+nunktMech;
        end
    end
end

for i=1:nface
    for j=1:order-1
        if unkfatp3(i,j)<0
            unkfatp3(i,j)=abs(unkfatp3(i,j))+nunkt+nunktMech;
        end
    end
end

sol(1+nunkt+nunktMech:npec+nunkt+nunktMech)=known(1:npec);
npecMech=unknown.Mech.npec;

nunktotal=nunkt+npec+nunktMech+npecMech;
% correct for zero terms fix all zero values to be at nunktotal+1

npecMech=unknown.Mech.npec;
nunktMech=unknown.Mech.nunkt;
nunktSystem=unknown.system.nunkt;
npecSystem=npec+npecMech;
%helpSol(1:nunktMech+npecMech)=sol(nunktotal+1:end);
%helpSol=helpSol';
%sol(nunkt+npec+1)=zeros(1,nrhs);

for i=1:nside
    if unkz(i)==0
        unkz(i)=nunktotal+1;
    end
end

for i=1:nside
    for j=1:order
        if unksid(i,j)==0
            unksid(i,j)=nunktotal+1;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if unkfatp1(i,j)==0
            unkfatp1(i,j)=nunktotal+1;
        end
    end
end

for i=1:nface
    for j=1:(order*order-order)/2
        if(unkfatp2(i,j)==0)
            unkfatp2(i,j)=nunktotal+1;
        end
    end
end

for i=1:nface
    for j=1:order-1
        if unkfatp3(i,j)==0
            unkfatp3(i,j)=nunktotal+1;
        end
    end
end

if order>=3
for i=1:nelem
    for j=1:(order-2)*(order-1)*(order+1)/2
        if unkint(i,j)==0
            unkint(i,j)=nunktotal+1;
        end
    end
end
end

% Save relevant data to unknown structure
unknown.EM.unkz=unkz;
unknown.EM.unksid=unksid;
unknown.EM.unkint=unkint;
unknown.EM.unkfatp1=unkfatp1;
unknown.EM.unkfatp2=unkfatp2;
unknown.EM.unkfatp3=unkfatp3;

%========================================================================================
% Include Dirichlet values for Mechanical problem
%========================================================================================

% Extract relevant data from unknown structure

unkvertexMech=unknown.Mech.unkvertexMech;
unkvertexSystem=unknown.system.unkvertexSystem;
unkedgesMech=unknown.Mech.unkedgesMech;
unkedgesSystem=unknown.system.unkedgesSystem;
unkfacesMech=unknown.Mech.unkfacesMech;
unkfacesSystem=unknown.system.unkfacesSystem;
unkinteriorsMech=unknown.Mech.unkinteriorsMech;
unkinteriorsSystem=unknown.system.unkinteriorsSystem;

 nelem=mesh.mech.Nelements;
nside=mesh.mech.edge.nside;
nface=mesh.mech.face.nface;
nNodes=mesh.mech.nNodes;







% Vertex unknowns

for i=1:nNodes
    for k=1:3
        if unkvertexMech(i,k)<0
            unkvertexMech(i,k)=abs(unkvertexMech(i,k))+nunktMech;
            unkvertexSystem(i,k)=abs(unkvertexSystem(i,k))+nunktSystem;
        end
    end
end

% Edges unknowns
for j=0:orderH1-2
    for i=1:nside
        for k=1:3
            if unkedgesMech(i,j+1,k)<0
                unkedgesMech(i,j+1,k)=abs(unkedgesMech(i,j+1,k))+nunktMech;
                unkedgesSystem(i,j+1,k)=abs(unkedgesSystem(i,j+1,k))+nunktSystem;
            end
        end
    end
end

% Face unknowns

   for i =1:nface
            for j=1:((orderH1-1)^2-(orderH1-1))/2
                for k=1:3
                    if unkfacesMech(i,j,k)<0
                        unkfacesMech(i,j,k)=abs(unkfacesMech(i,j,k))+nunktMech;
                        unkfacesSystem(i,j,k)=abs(unkfacesSystem(i,j,k))+nunktSystem;
                    end
                end
            end
   end
   
% Interior unknowns
   
      for i = 1:nelem
       
            k=0;
            for ii=0:orderH1-4
                for jj=0:orderH1-4
                    for kk=0:orderH1-4
                        if ii+jj+kk <= orderH1-4
                            k=k+1;
                            for l=1:3
                                if unkinteriorsMech(i,k,l)<0
                                    unkinteriorsMech(i,k,l)=abs(unkinteriorsMech(i,k,l))+nunktMech;
                                    unkinteriorsSystem(i,k,l)=abs(unkinteriorsSystem(i,k,l))+nunktSystem;
                                end
                            end
                        end
                    end
                end
            end
      end
                
nunkTotal=length(sol);

%sol(nunkt+npec+1:nunkt+npec+nunktMech+npecMech)=helpSol(1:end);
sol(nunkt+npec+nunktMech+1:nunkt+npec+nunktMech+npecMech)=known_mech(1:end);
%nunkFinal=length(sol);
%sol(nunkFinal+1)=0;


% All zero values at final position


% Vertex unknowns

for i=1:nNodes
    for k=1:3
        if unkvertexMech(i,k)==0
            unkvertexMech(i,k)=nunktotal+1;
            unkvertexSystem(i,k)=nunktotal+1;
        end
    end
end

% Edges unknowns
for j=0:orderH1-2
    for i=1:nside
        for k=1:3
            if unkedgesMech(i,j+1,k)==0
                unkedgesMech(i,j+1,k)=nunktotal+1;
                unkedgesSystem(i,j+1,k)=nunktotal+1;
            end
        end
    end
end

% Face unknowns

   for i =1:nface
            for j=1:((orderH1-1)^2-(orderH1-1))/2
                for k=1:3
                    if unkfacesMech(i,j,k)==0
                        unkfacesMech(i,j,k)=nunktotal+1;
                        unkfacesSystem(i,j,k)=nunktotal+1;
                    end
                end
            end
   end
   
% Interior unknowns
   
      for i = 1:nelem
       
            k=0;
            for ii=0:orderH1-4
                for jj=0:orderH1-4
                    for kk=0:orderH1-4
                        if ii+jj+kk <= orderH1-4
                            k=k+1;
                            for l=1:3
                                if unkinteriorsMech(i,k,l)==0
                                    unkinteriorsMech(i,k,l)=nunktotal+1;
                                    unkinteriorsSystem(i,k,l)=nunktotal+1;
                                end
                            end
                        end
                    end
                end
            end
      end
sol(nunktotal+1)=0;

unknown.Mech.unkvertexMech=unkvertexMech;
unknown.Mech.unkedgesMech=unkedgesMech;
unknown.Mech.unkfacesMech=unkfacesMech;
unknown.Mech.unkinteriorsMech=unkinteriorsMech;

unknown.system.unkvertexSystem=unkvertexSystem;
unknown.system.unkedgesSystem=unkedgesSystem;
unknown.system.unkfacesSystem=unkfacesSystem;
unknown.system.unkinteriorsSystem=unkinteriorsSystem;

















