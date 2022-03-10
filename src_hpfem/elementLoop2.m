function [mass,damp,damppre,stiff,Resid]=elementLoop2(Static,mesh,Basis,Quadrature,unknown,probdata,X,nUnk,nunkMech,gorder,FEspacesInfo,probstatic,known,known_mech,DirichletFaceBasis)
%=========================================================================
% Extract the data from the data structure
%=========================================================================
% Extract the relevant data from the mesh structure
intma = mesh.intma;
mat = mesh.mat;
nelem = mesh.Nelements;
edgecof = mesh.edgecof;
facecof = mesh.facecof;
coord = mesh.Coordinates;
glob = mesh.edge.glob;
globfa = mesh.face.globfa;
matc=mesh.matc;
subFlag=mesh.subFlag;
cond=mesh.face.cond;

% Extract relevant data from unknown structure
unkz=unknown.EM.unkz;
unksid=unknown.EM.unksid;
unkint=unknown.EM.unkint;
unkfatp1=unknown.EM.unkfatp1;
unkfatp2=unknown.EM.unkfatp2;

unkfatp3=unknown.EM.unkfatp3;
nunkt=unknown.EM.nunkt;
nunktMech=unknown.Mech.nunkt;

% Extract relevant data from probdata structure
order=probdata.jb.order;
sigma=probdata.matr.sigma;
omega=probdata.matr.omega;
muz=probdata.matr.muz;
mu=probdata.matr.mu;
lambda=probdata.matr.lambda;
G=probdata.matr.G;


regTerm=probdata.matr.regTerm;

% Extract relevant data from FEspacesInfo structure
orderH1=FEspacesInfo.orderH1;
esizet=FEspacesInfo.esizet;
esizeH1=FEspacesInfo.esizeH1;

% Extract relevant data from Quadrature structure

intfw=Quadrature.intfw;
intfxi=Quadrature.intfxi;
intfet=Quadrature.intfet;
nipf=Quadrature.nipf;

% Extract relevant data from DirichletFaceBasis structure

ephdfx1=DirichletFaceBasis.ephdfx1;
ephdfy1=DirichletFaceBasis.ephdfy1;
ephdfz1=DirichletFaceBasis.ephdfz1;
ephdfx2=DirichletFaceBasis.ephdfx2;
ephdfy2=DirichletFaceBasis.ephdfy2;
ephdfz2=DirichletFaceBasis.ephdfz2;

phh1f1=DirichletFaceBasis.phh1f1;
phh1f2=DirichletFaceBasis.phh1f2;
phH1f1=DirichletFaceBasis.phH1f1;
phH1f2=DirichletFaceBasis.phH1f2;

gphfx1=DirichletFaceBasis.gphfx1;
gphfy1=DirichletFaceBasis.gphfy1;
gphfz1=DirichletFaceBasis.gphfz1;
gphfx2=DirichletFaceBasis.gphfx2;
gphfy2=DirichletFaceBasis.gphfy2;
gphfz2=DirichletFaceBasis.gphfz2;



%=========================================================================
% Begin initialisations
%=========================================================================
% Define the size of the matrix
sizeMat    = unknown.EM.nunkt;


% Initialise arrays ready for assembly
Resid      = zeros(sizeMat+nunktMech,1);

 I          = zeros(nUnk+nunkMech,1);
 J          = zeros(nUnk+nunkMech,1);


Kval       = zeros(nUnk+nunkMech,1);
Cval       = zeros(nUnk+nunkMech,1);
Cvalpre    = zeros(nUnk+nunkMech,1);
Mval       = zeros(nUnk+nunkMech,1);


% Initialise variables
count      = 0;
nmst =0;


Volume=0;
AreaInt=0;
%=========================================================================
% Loop through elements
%=========================================================================
for i=1:nelem
    mue=mu(mat(i));
    if matc(i)==1
        if probstatic==0
        %     conductor
        kappalow=sigma(mat(i))*muz;
        elseif probstatic==1
            kappalow=-complex(0,1)*regTerm*muz;
            
        end
    else
       
        %     regularised free space
        kappalow=-complex(0,1)*regTerm*muz;
        
       
    end
    
    if probdata.sol.regopt ==1
        if matc(i)~=1
            kappa=0;
        else
            kappa=kappalow;
        end
    else
        kappa=kappalow;
    end
    
    % Obtain the global DOF numbers for this element
        bhelp=rowfun(unknown,order,orderH1,mesh,i,esizet,esizeH1,subFlag(i));
    
       % transfer coordinates to local array
    xy = coord(intma(i,1:4),1:3);
    x(i)=(coord(intma(i,1),1)+coord(intma(i,2),1)+coord(intma(i,3),1)+coord(intma(i,4),1))/4;
    y(i)=(coord(intma(i,1),2)+coord(intma(i,2),2)+coord(intma(i,3),2)+coord(intma(i,4),2))/4;
    z(i)=(coord(intma(i,1),3)+coord(intma(i,2),3)+coord(intma(i,3),3)+coord(intma(i,4),3))/4;
    flag=0 ;
    gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
    mycoord = zeros(gesizet,3);
    %transfer coefficents to locations vertices
    mycoord(1:4,1:3) = xy(1:4,1:3);
    lec=zeros(6,gorder,3);
    lfc=zeros(4,gorder*(gorder-1)/2,3);
    if gorder > 0
        for j=1:6
            for p=1:gorder
                for k=1:3
                    lec(j,p,k)=edgecof(glob(i,j),((p-1)*3)+k);
                    if abs(edgecof(glob(i,j),((p-1)*3)+k))>0.00000001
                        flag=1;
                    end
                    mycoord(4+j+6*(p-1),k)=lec(j,p,k);
                end
            end
        end
        
        for j=1:4
            for p=1:gorder*(gorder-1)/2
                for k=1:3
                    lfc(j,p,k)=facecof(globfa(i,j),((p-1)*3)+k);
                    mycoord(4+6*gorder+(j-1)*gorder*(gorder-1)/2+p,k)= lfc(j,p,k);
                end
            end
        end
    end
    
    % Display the number of elements completed at intervals in the run
    count=count+1;
    if count==1000
        count=0;
        disp(['number of elements complete, nelem = ',num2str(i)]);
    end
    
 
    
    eltype=mesh.eltype(i);
    material=mat(i);
    
    %---------------------------------------------------------------------
    % Determine the elemental solution vector
    %---------------------------------------------------------------------
    % Extract the static solution
    Adc=zeros(esizet,1);
      A=zeros(esizet,1);
    dAdt=zeros(esizet,1);
    dAdt2=zeros(esizet,1);
    u=zeros(3*esizeH1,1);
    udc=zeros(3*esizeH1,1);
    dudt=zeros(3*esizeH1,1);
    dudt2=zeros(3*esizeH1,1);
    nunktEM=unknown.EM.nunkt;
    nunkt=unknown.system.nunkt;
    npec=unknown.EM.npec;
    
   
      for ii=1:esizet
        row=bhelp(ii);
        if row>0
            Adc(ii,1)=Static.sol(row,1);
               A(ii,1)=X(row,1);
            dAdt(ii,1)=X(row,2);
            dAdt2(ii,1)=X(row,3);
        elseif row<0
              Adc(ii,1)=Static.sol(abs(row)+nunktEM,1);
               A(ii,1)=X(abs(row)+nunktEM,1);
            dAdt(ii,1)=X(abs(row)+nunktEM,2);
            dAdt2(ii,1)=X(abs(row)+nunktEM,3);
        end
            
      end
      
       for ii=esizet+1:esizet+3*esizeH1
        row=bhelp(ii);
        if row>0
            udc(ii-esizet,1)=Static.sol(row+npec,1);
               u(ii-esizet,1)=X(row+npec,1);
            dudt(ii-esizet,1)=X(row+npec,2);
            dudt2(ii-esizet,1)=X(row+npec,3);
        elseif row<0
              udc(ii-esizet,1)=Static.sol(abs(row)+nunkt,1);
               u(ii-esizet,1)=X(abs(row)+nunkt,1);
            dudt(ii-esizet,1)=X(abs(row)+nunkt,2);
            dudt2(ii-esizet,1)=X(abs(row)+nunkt,3);
        end
      end
    

    % Store the dynamic solution in the array
    Xa    = [A dAdt dAdt2];
    Xu    = [u dudt dudt2];
    
    %---------------------------------------------------------------------
    % Determine the Elemental matrices
    %---------------------------------------------------------------------
    % Compute the elemental matrices by gauss integration
    
    [M,C,K,Res,Volume,AreaInt]=GaussInteg2(Basis,DirichletFaceBasis,Quadrature,probdata,flag,gorder,mue,esizet,esizeH1,omega,Adc,Xa,Xu,mycoord,eltype,probstatic,material,subFlag(i),i,lec,lfc,xy,cond,matc,Volume,AreaInt);
    % Regularization
    
      C1(1:6,1:6)=kappalow*C(1:6,1:6);
    Cpre(1:6,1:6)=C1(1:6,1:6);
    
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
nintbas=k;
 
    
    %     higher order block
    for p=1:esizet
        if p<= 6
            C1(p,7:esizet)=kappa*C(p,7:esizet);
        else
            C1(p,1:esizet)=kappa*C(p,1:esizet);
            Cpre(p,1:esizet)=abs(kappa)*C(p,1:esizet);
        end
    end
    
    if probdata.sol.regopt==1
        % high order edges (gradients)
        Cpre(7:6*(order+1),7:6*(order+1))=-complex(0,1)*abs(kappa)*C(7:6*(order+1),7:6*(order+1));
        
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappa)*C(st:en,st:en);
        
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappalow)*C(st:en,st:en);
        
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappa)*C(st:en,st:en);
       
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappalow)*C(st:en,st:en);
    else
        % high order edges (gradients)
        Cpre(7:6*(order+1),7:6*(order+1))=-complex(0,1)*abs(kappalow)*C(7:6*(order+1),7:6*(order+1));
        
        % high order faces (gradients)
        st=6*(order+1)+1;en=6*(order+1)+4*(order*order-order)/2;
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappalow)*C(st:en,st:en);
        
        % high order faces (non-gradients)
        st=6*(order+1)+4*(order*order-order)/2+1; en=6*(order+1)+4*(2*(order*order-order)/2+(order-1));
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappalow)*C(st:en,st:en);
        
        % high order interiors (gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+1;
        en=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas;
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappalow)*C(st:en,st:en);
        
        % high order interiors (non-gradients)
        st=6*(order+1)+4*(2*(order*order-order)/2+(order-1))+nintbas+1;
        en=esizet;
        Cpre(st:en,st:en)=-complex(0,1)*abs(kappalow)*C(st:en,st:en);
    end
    
        if eltype==1
        phH1f=phH1f1;
        phh1f=phh1f1;
        gphfx=gphfx1;
        gphfy=gphfy1;
        gphfz=gphfz1;
        ephdfx=ephdfx1;
        ephdfy=ephdfy1;
        ephdfz=ephdfz1;
        elseif eltype==2
        phH1f=phH1f2;
        phh1f=phh1f2;
        gphfx=gphfx2;
        gphfy=gphfy2;
        gphfz=gphfz2;
        ephdfx=ephdfx2;
        ephdfy=ephdfy2;
        ephdfz=ephdfz2;
        end
        

        

    
%     Res= NeumannBC(Res,xy,esizet,esizeH1,intfxi,intfet,intfw,nipf,eltype,cond,...
%     i,lec,lfc,flag,gorder,probdata,...
%     mycoord,phh1f,phH1f,gphfx,gphfy,gphfz,ephdfx,ephdfy,ephdfz,mue,probstatic,subFlag(i),material);
%     

C1(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1)=C(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1);
Cpre(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1)=C(esizet+1:esizet+3*esizeH1,esizet+1:esizet+3*esizeH1);
    

    %---------------------------------------------------------------------
    % Linear system assembly in vector format
    %---------------------------------------------------------------------
    [Mval,Cval,Cvalpre,Kval,Resid,I,J,nmst]=matrixAssembly(K,C1,Cpre,M,Res,Xa,Mval,Cval,Cvalpre,Kval,Resid,I,J,bhelp,probstatic,esizet,esizeH1,known,known_mech,nmst,omega);
    
end


%=========================================================================
% Construct the Linear system in sparse format
%=========================================================================
% Build stiffness and damping matrices from I,J,Kval and Cval
stiff = sparse(I(1:nmst),J(1:nmst),Kval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damp  = sparse(I(1:nmst),J(1:nmst),Cval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
mass  = sparse(I(1:nmst),J(1:nmst),Mval(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);
damppre=sparse(I(1:nmst),J(1:nmst),Cvalpre(1:nmst),nunktEM+nunktMech,nunktEM+nunktMech);




% Clear the system vectors
clear I
clear J
clear Kval
clear Cval
clear Mval


end