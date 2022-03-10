% Function for outputing the solution on sub-divided tetrahedrons for
% purpose of visualisation in Paraview.

% Two sets of files are created for the given "job" name specified in the
% problem file
% 1) .vtm block file specifying the details of the vtu files for each different subdomain
% 2) .vtu specifiying the detials of the solution and subdivided grid for
% a particular subdomain - one is created for each subdomain

function pointvalue_multi(mesh,unknown,probdata,Solution,probstatic,freq)

% Extract relevant data from mesh structure
glob=mesh.edge.glob;
globfa=mesh.face.globfa;
eltype=mesh.eltype;
coord=mesh.Coordinates;
intma=mesh.intma;
mat=mesh.mat;
edgecof=mesh.edgecof;
facecof=mesh.facecof;
subFlag=mesh.subFlag;
sol=Solution;

% Extract relevant data from FEspacesInfo structure
order=probdata.order;
esizet=probdata.esizet;
orderH1=probdata.orderH1;
esizeH1=probdata.esizeH1;
gorder=probdata.jb.gorder;
job=probdata.jb.job;




% split the mesh and write out as multi-blocks

omega=2*pi*freq;

[nelem ~]=size(intma);

phxstore1 = [];
phystore1 = [];
phzstore1 = [];
phxstore2 = [];
phystore2 = [];
phzstore2 = [];

cphxstore1 = [];
cphystore1 = [];
cphzstore1 = [];
cphxstore2 = [];
cphystore2 = [];
cphzstore2 = [];

gphxstore1=[];
gphystore1=[];
gphzstore1=[];
gphxstore2=[];
gphystore2=[];
gphzstore2=[];
gphxstoreQuad=[];
gphystoreQuad=[];
gphzstoreQuad=[];

phh1store1=[];
phh1store2=[];
H1basisStore1=[];
H1basisStore2=[];

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

% Split the tetrahedron:
% produce a splitting of the reference element
% problems with MATLAB's delauney so fix number of points to be 4
% 3 does not produce correct triangulation!
if order <=4%5%2
   [Splitref,intmasplit] = nodalCoordinatesTetraMinus1(order+1);
   %[Splitref,intmasplit] = nodalCoordinatesTetraMinus1(1);
else
    [Splitref,intmasplit] = nodalCoordinatesTetraMinus1(4);
end
[netet dum]=size(intmasplit);

[nppt dum]=size(Splitref);

gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
for i = 1:nppt
    xi = Splitref(i,1);
    eta = Splitref(i,2);
    zeta = Splitref(i,3);
    ph =  basis(xi,eta,zeta,axi,aeta,azeta,order,1,esizet);
    phxstore1 = [phxstore1,ph(:,1)];
    phystore1 = [phystore1,ph(:,2)];
    phzstore1 = [phzstore1,ph(:,3)];
    ph =  basis(xi,eta,zeta,axi,aeta,azeta,order,2,esizet);
    phxstore2 = [phxstore2,ph(:,1)];
    phystore2 = [phystore2,ph(:,2)];
    phzstore2 = [phzstore2,ph(:,3)];
    ph =  curlbasis(xi,eta,zeta,asxi,aseta,aszeta,order,1,esizet);
    cphxstore1 = [cphxstore1,ph(:,1)];
    cphystore1 = [cphystore1,ph(:,2)];
    cphzstore1 = [cphzstore1,ph(:,3)];
    ph =  curlbasis(xi,eta,zeta,asxi,aseta,aszeta,order,2,esizet);
    cphxstore2 = [cphxstore2,ph(:,1)];
    cphystore2 = [cphystore2,ph(:,2)];
    cphzstore2 = [cphzstore2,ph(:,3)];
    % geometry
    gph=gbasish1(xi,eta,zeta,axi,aeta,azeta,gorder+1,1,gesizet);
    gphxstore1=[gphxstore1, gph(:,1)];
    gphystore1=[gphystore1, gph(:,2)];
    gphzstore1=[gphzstore1, gph(:,3)];
    gph=gbasish1(xi,eta,zeta,axi,aeta,azeta,gorder+1,2,gesizet);
    gphxstore2=[gphxstore2, gph(:,1)];
    gphystore2=[gphystore2, gph(:,2)];
    gphzstore2=[gphzstore2, gph(:,3)];
    gphQuad=BasisQuadraticGeom(xi,eta,zeta);
    gphxstoreQuad=[gphxstoreQuad, gphQuad(:,1)];
    gphystoreQuad=[gphystoreQuad, gphQuad(:,2)];
    gphzstoreQuad=[gphzstoreQuad, gphQuad(:,3)];
    phh1=basish1(gesizet,xi,eta,zeta,gorder+1,1);
    phh1store1=[phh1store1, phh1(:,1)];
    phh1=basish1(gesizet,xi,eta,zeta,gorder+1,2);
    phh1store2=[phh1store2, phh1(:,1)];
    th1=basish1(esizeH1,xi,eta,zeta,orderH1,1);
    H1basisStore1=[H1basisStore1,th1(:,1)];
    
    th1=basish1(esizeH1,xi,eta,zeta,orderH1,2);
    H1basisStore2=[H1basisStore2,th1(:,1)];
    
end

% for each sub-domain find the list of elements contained
nmat=max(mat);
helpmat=zeros(nmat,nelem);
helpmate=zeros(nmat);
for i=1:nelem
    helpmate(mat(i))=helpmate(mat(i))+1;
    helpmat(mat(i),helpmate(mat(i)))=i;
end
% if gorder==0 we need to define dummy arrays:
lec=[];
lfc=[];
for imat = 1:nmat
    disp(['Processing material ',num2str(imat)]);
    nelemn=netet*helpmate(imat);
    npoinn=0;
    intman=[];
    coordn=zeros(nelem*nppt,3);
    pen=zeros(nelem*nppt,3);
    phen=pen;
    displ=zeros(nelem*nppt,3);
    displx=zeros(nelem*nppt,3);
    disply=zeros(nelem*nppt,3);
    displz=zeros(nelem*nppt,3);
    for je=1:helpmate(imat)
        i=helpmat(imat,je);
        xy = coord(intma(i,1:4),1:3);
        
    gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
    
    mycoord = zeros(gesizet,3);
    
    %transfer coefficents to locations vertices
    mycoord(1:4,1:3) = xy(1:4,1:3);
    flag=0;
    
    if gorder>1
    for j=1:6
        for p=1:gorder
            for k=1:3
                lec(j,p,k)=edgecof(glob(i,j),((p-1)*3)+k);
                if abs(edgecof(glob(i,j),((p-1)*3)+k))>0.0000001
                    flag=1;
                end
            end
        end
    end
    
    for j=1:4
        for p=1:gorder*(gorder-1)/2
            for k=1:3
                lfc(j,p,k)=facecof(globfa(i,j),((p-1)*3)+k);
            end
        end
    end

    
    % edge functions
    for ii=1:6
        for p=1:gorder
            for j=1:3
                
                mycoord(4+ii+6*(p-1),j)=lec(ii,p,j);
            end
        end
    end
    
    % face functions
    for iii=1:4
        for ii=1:(gorder-1)*gorder/2
            for j=1:3
                mycoord(4+6*gorder+(iii-1)*gorder*(gorder-1)/2+ii,j)= lfc(iii,ii,j);
            end
        end
    end
    elseif gorder==1
        mycoord=mesh.mycoord(:,:,i);
        flag=1;
    end
    
        
        %         bhelp=rowfun(unknown,order,orderH1,mesh,...
        %             i,esizet,esizeH1,subFlag(i));
        
        bhelp=unknown.EM.unknowns(:,i);
        bhelp2=unknown.Mech.unknownsD(:,i);
        % generate connectivities
        intman( (je-1)*netet+1:je*netet,:)= npoinn+intmasplit;
        
        
        % Stored basis functions
        if eltype(i)==1
            cphxstore=cphxstore1;
            cphystore=cphystore1;
            cphzstore=cphzstore1;
            phxstore=phxstore1;
            phystore=phystore1;
            phzstore=phzstore1;
            gphxstore=gphxstore1;
            gphystore=gphystore1;
            gphzstore=gphzstore1;
            phh1store=phh1store1;
            H1basisStore=H1basisStore1;
        else
            cphxstore=cphxstore2;
            cphystore=cphystore2;
            cphzstore=cphzstore2;
            phxstore=phxstore2;
            phystore=phystore2;
            phzstore=phzstore2;
            gphxstore=gphxstore2;
            gphystore=gphystore2;
            gphzstore=gphzstore2;
            phh1store=phh1store2;
            H1basisStore=H1basisStore2;
        end
        he1=zeros(esizet,1);
        he1(bhelp>0) = sol(bhelp(bhelp>0),1);
        
        dispCoeff=zeros(3*esizeH1,1);
        if subFlag(i)==2
         dispCoeff(bhelp2>0)=sol(bhelp2(bhelp2>0),1);   
        end
        
        if flag==0
            % straight sided element
            % constant Jacobian
            gph=[gphxstore(:,1) gphystore(:,1) gphzstore(:,1)];
            [axi,aeta,azeta,asxi,aseta,aszeta,~]=jacobian_pre(flag,gesizet,gph,mycoord);
            
            
            for j = 1:nppt
                ph1=phh1store(:,j);
                % computx x,y,z
                [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
                
                npoinn=npoinn+1;
                coordn(npoinn,:) = [x,y,z];
                
                phh=(cphxstore(:,j)*asxi(1:3))+...
                    (cphystore(:,j)*aseta(1:3))+...
                    (cphzstore(:,j)*aszeta(1:3));
                
                
                hen1 = phh'*he1/probdata.matr.muz/probdata.matr.mu(mat(i));
                
                %       store the magentic field at each point
                phen(npoinn,:) = hen1.';
                
%                 %==========================================================
%                 % Store the magnetic field coming from Biot-Savart law
%                 %==========================================================
%                 % ---------------------------------------------------------
%                 % Compute B for Coil 1 (right)
%                 %----------------------------------------------------------
%                 zCenter=1.6485e-1;
%                 z1=(194.7e-3)-zCenter;
%                 zB=z-zCenter;
%                 rho1=250e-3;
%                 rho2=288.4e-3;
%                 rho=sqrt(x^2+y^2);
%                 JDC=253.2*1e6*1.256637061435917e-06;
%                 B1(npoinn,1)=0.5*JDC*((-rho)*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1))+sqrt((rho-rho2)^2+(zB-z1)^2)-sqrt((rho-rho1)^2+(zB-z1)^2)-sqrt((rho-rho2)^2+(zB+z1)^2)+sqrt((rho-rho1)^2+(zB+z1)^2));
%                 B1(npoinn,3)=-0.5*JDC*(-rho*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+zB+z1)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+zB+z1))+(zB-z1)*log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-(zB-z1)*log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-(zB+z1)*log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1)+(zB+z1)*log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2));
%                 phi=atan2(y,x);
%                 B3(npoinn,1)=cos(phi)*B1(npoinn,1);
%                 B3(npoinn,2)=sin(phi)*B1(npoinn,1);
%                 B3(npoinn,3)=B1(npoinn,3);
%                 %----------------------------------------------------------
%                 % Compute B for Coil 2 (left)
%                 %----------------------------------------------------------
%                 zCenter=1.6485e-1;
%                 z1=(194.7e-3)-zCenter;
%                 zB=z+zCenter;
%                 rho1=250e-3;
%                 rho2=288.4e-3;
%                 rho=sqrt(x^2+y^2);
%                 JDC=253.2*1e6*1.256637061435917e-06;
%                 B2(npoinn,1)=0.5*JDC*((-rho)*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1))+sqrt((rho-rho2)^2+(zB-z1)^2)-sqrt((rho-rho1)^2+(zB-z1)^2)-sqrt((rho-rho2)^2+(zB+z1)^2)+sqrt((rho-rho1)^2+(zB+z1)^2));
%                 B2(npoinn,3)=-0.5*JDC*(-rho*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+zB+z1)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+zB+z1))+(zB-z1)*log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-(zB-z1)*log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-(zB+z1)*log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1)+(zB+z1)*log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2));
%                 phi=atan2(y,x);
%                 B4(npoinn,1)=cos(phi)*B2(npoinn,1);
%                 B4(npoinn,2)=sin(phi)*B2(npoinn,1);
%                 B4(npoinn,3)=B2(npoinn,3);
%                 B1(npoinn,:)=B3(npoinn,:)+B4(npoinn,:);
%                 %==========================================================
                
                phh = (phxstore(:,j)*axi(1:3))+...
                    (phystore(:,j)*aeta(1:3))+...
                    (phzstore(:,j)*azeta(1:3));
                
                en1 =complex(0,-1)*omega*phh'*he1*probdata.matr.sigma(mat(i));
                
                %       store the eddy current at each point
                pen(npoinn,:) =en1.';
                if subFlag(i)==2
                    H1bas(1:esizeH1)=H1basisStore(1:esizeH1,j);
                    
                    H1bas3D=zeros(3,3*esizeH1);
                    H1bas3D(1,1:3:end)=H1bas;
                    H1bas3D(2,2:3:end)=H1bas;
                    H1bas3D(3,3:3:end)=H1bas;
                    
                    
                    dispPoint=H1bas3D*dispCoeff;
                    displ(npoinn,:)=dispPoint.';
                    displx(npoinn,1)=displ(npoinn,1);
                    displx(npoinn,2)=0;
                    displx(npoinn,3)=0;
                    disply(npoinn,1)=0;
                    disply(npoinn,2)=displ(npoinn,2);
                    disply(npoinn,3)=0;
                    displz(npoinn,1)=0;
                    displz(npoinn,2)=0;
                    displz(npoinn,3)=displ(npoinn,3);
                else
                    displ(npoinn,:)=zeros(1,3);
                    displx(npoinn,1)=0;
                    displx(npoinn,2)=0;
                    displx(npoinn,3)=0;
                    disply(npoinn,1)=0;
                    disply(npoinn,2)=0;
                    disply(npoinn,3)=0;
                    displz(npoinn,1)=0;
                    displz(npoinn,2)=0;
                    displz(npoinn,3)=0;
                end
                
            end
        else
            % curved sided element
            
            
            for j = 1:nppt
                ph1=phh1store(:,j);
                    xi = Splitref(j,1);
                    eta = Splitref(j,2);
                    zeta = Splitref(j,3);
                % computx x,y,z
                if gorder>1
                    [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
                else
                    [x,y,z]= getxyzq(mycoord,xi,eta,zeta);
                end
                npoinn=npoinn+1;
                coordn(npoinn,:) = [x,y,z];
                
                gph=[gphxstore(:,j) gphystore(:,j) gphzstore(:,j)];
                
                   if gorder>1 
                    [axi,aeta,azeta,asxi,aseta,aszeta,~]=jacobian_pre(flag,gesizet,gph,mycoord) ;
                   else
                    gphQuad=[gphxstoreQuad(:,j) gphystoreQuad(:,j) gphzstoreQuad(:,j)];
                                
                    % evaluate covairant mapping
                    [axi,aeta,azeta,asxi,aseta,aszeta,~]=jacobian_pre(flag,gesizet,gphQuad,mycoord);
                   end
                
                phh=(cphxstore(:,j)*asxi(1:3))+...
                    (cphystore(:,j)*aseta(1:3))+...
                    (cphzstore(:,j)*aszeta(1:3));
                
                
                hen1 = phh'*he1/probdata.matr.muz/probdata.matr.mu(mat(i));
                
                
                %       store the magentic field at each point
                phen(npoinn,:) = hen1.';
                
%                 %=========================================================
%                 % Store the magnetic field coming from Biot-Savart law
%                 %=========================================================
%                 % ---------------------------------------------------------
%                 % Compute B for Coil 1 (right)
%                 %----------------------------------------------------------
%                 zCenter=1.6485e-1;
%                 z1=(194.7e-3)-zCenter;
%                 zB=z-zCenter;
%                 rho1=250e-3;
%                 rho2=288.4e-3;
%                 rho=sqrt(x^2+y^2);
%                 JDC=253.2*1e6*1.256637061435917e-06;
%                 B1(npoinn,1)=0.5*JDC*((-rho)*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1))+sqrt((rho-rho2)^2+(zB-z1)^2)-sqrt((rho-rho1)^2+(zB-z1)^2)-sqrt((rho-rho2)^2+(zB+z1)^2)+sqrt((rho-rho1)^2+(zB+z1)^2));
%                 B1(npoinn,3)=-0.5*JDC*(-rho*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+zB+z1)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+zB+z1))+(zB-z1)*log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-(zB-z1)*log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-(zB+z1)*log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1)+(zB+z1)*log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2));
%                 phi=atan2(y,x);
%                 B3(npoinn,1)=cos(phi)*B1(npoinn,1);
%                 B3(npoinn,2)=sin(phi)*B1(npoinn,1);
%                 B3(npoinn,3)=B1(npoinn,3);
%                 %----------------------------------------------------------
%                 % Compute B for Coil 2 (left)
%                 %----------------------------------------------------------
%                 zCenter=1.6485e-1;
%                 z1=(194.7e-3)-zCenter;
%                 zB=z+zCenter;
%                 rho1=250e-3;
%                 rho2=288.4e-3;
%                 rho=sqrt(x^2+y^2);
%                 JDC=253.2*1e6*1.256637061435917e-06;
%                 B2(npoinn,1)=0.5*JDC*((-rho)*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1))+sqrt((rho-rho2)^2+(zB-z1)^2)-sqrt((rho-rho1)^2+(zB-z1)^2)-sqrt((rho-rho2)^2+(zB+z1)^2)+sqrt((rho-rho1)^2+(zB+z1)^2));
%                 B2(npoinn,3)=-0.5*JDC*(-rho*(log(sqrt((rho-rho2)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho1)^2+(zB-z1)^2)+zB-z1)-log(sqrt((rho-rho2)^2+(zB+z1)^2)+zB+z1)+log(sqrt((rho-rho1)^2+(zB+z1)^2)+zB+z1))+(zB-z1)*log(sqrt((rho-rho1)^2+(zB-z1)^2)+rho-rho1)-(zB-z1)*log(sqrt((rho-rho2)^2+(zB-z1)^2)+rho-rho2)-(zB+z1)*log(sqrt((rho-rho1)^2+(zB+z1)^2)+rho-rho1)+(zB+z1)*log(sqrt((rho-rho2)^2+(zB+z1)^2)+rho-rho2));
%                 phi=atan2(y,x);
%                 B4(npoinn,1)=cos(phi)*B2(npoinn,1);
%                 B4(npoinn,2)=sin(phi)*B2(npoinn,1);
%                 B4(npoinn,3)=B2(npoinn,3);
%                 B1(npoinn,:)=B3(npoinn,:)+B4(npoinn,:);
%                 %==========================================================
                
                
                phh = (phxstore(:,j)*axi(1:3))+...
                    (phystore(:,j)*aeta(1:3))+...
                    (phzstore(:,j)*azeta(1:3));
                
                
                en1 =complex(0,-1)*omega*phh'*he1*probdata.matr.sigma(mat(i));
                
                %       store the eddy current at each point
                pen(npoinn,:) =en1.';
                
                if subFlag(i)==2
                    H1bas(1:esizeH1)=H1basisStore(1:esizeH1,j);
                    
                    H1bas3D=zeros(3,3*esizeH1);
                    H1bas3D(1,1:3:end)=H1bas;
                    H1bas3D(2,2:3:end)=H1bas;
                    H1bas3D(3,3:3:end)=H1bas;
                    
                    dispPoint=H1bas3D*dispCoeff;
                    displ(npoinn,:)=dispPoint.';
                    displx(npoinn,1)=displ(npoinn,1);
                    displx(npoinn,2)=0;
                    displx(npoinn,3)=0;
                    disply(npoinn,1)=0;
                    disply(npoinn,2)=displ(npoinn,2);
                    disply(npoinn,3)=0;
                    displz(npoinn,1)=0;
                    displz(npoinn,2)=0;
                    displz(npoinn,3)=displ(npoinn,3);
                else
                    displ(npoinn,:)=zeros(1,3);
                    displx(npoinn,1)=0;
                    displx(npoinn,2)=0;
                    displx(npoinn,3)=0;
                    disply(npoinn,1)=0;
                    disply(npoinn,2)=0;
                    disply(npoinn,3)=0;
                    displz(npoinn,1)=0;
                    displz(npoinn,2)=0;
                    displz(npoinn,3)=0;
                    
                end
                
                
            end
            
        end
        
    end
    % write out this VTU file
    disp('writing data to the file...')
    if probstatic==0
        filename=[job '_SeveralIter' num2str(imat-1) '_' num2str(freq) 'Hz' '_p' num2str(orderH1) '.vtu']
        vtuk_puvw_write ( filename, npoinn, nelemn, ...
            coordn, intman, pen, phen,displ,displx,disply,displz)
        disp('done')
    else
        filename=[job '_StaticFine' num2str(imat-1) '_p' num2str(orderH1) '.vtu']
        vtuk_puvw_write ( filename, npoinn, nelemn, ...
            coordn, intman, pen,phen,displ)
        disp('done')
    end
end
% write out multi-block file
filename=[job '.vtm'];
output_unit = fopen(filename, 'w');

fprintf ( output_unit, '<VTKFile type="vtkMultiBlockDataSet" version="0.1" byte_order="BigEndian"\n');
fprintf ( output_unit, 'compressor="vtkZLibDataCompressor">\n');
fprintf ( output_unit, '<vtkMultiBlockDataSet>\n');
for imat=1:nmat
    filename=[job '_' num2str(imat-1) '.vtu'];
    text=['<DataSet group="' num2str(imat-1) '" dataset="' num2str(imat-1) '" file="' filename '"/>\n'];
    fprintf ( output_unit, text);
    
end
fprintf ( output_unit, '</vtkMultiBlockDataSet>\n');
fprintf ( output_unit, '</VTKFile>\n');
fclose(output_unit);

