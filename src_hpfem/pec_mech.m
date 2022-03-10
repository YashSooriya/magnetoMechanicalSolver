function [dirval,DirichletFaceBasis]= DirichletMech(mesh,unknown,Quadrature,gorder,order,orderH1,esizet,esizeH1,nipf,nrhs,probdata,probstatic)


% Extract relevant data form mesh structure
nelem=mesh.mech.Nelements;
nside=mesh.mech.nside;
coord=mesh.mech.Coordinates;
glob=mesh.mech.glob;
globfa=mesh.mech.globfa;
intma=mesh.mech.intma;
edgecof=mesh.mech.edgecof;
facecof=mesh.mech.facecof;
eltype=mesh.mech.eltype;
bcedge=mesh.mech.bcedge;
cond=mesh.mech.cond;
nNodes=mesh.mech.nNodes;
mapL2G_n=mesh.mech.mapL2G_n;

% Extract relevant data from unknown structure
unkvertexMech=unknown.unkvertexMech;
unkedgesMech=unknown.unkedgesMech;
npec=unknown.npec;
unkfacesMech=unknown.unkfacesMech;


% Extract relevant data from Quadrature structure
intfxi=Quadrature.intfxi;
intfet=Quadrature.intfet;
intfw=Quadrature.intfw;

probFlag=probdata.probFlag;

lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);
dirval = zeros(npec,1);

% find integration points along edge
% int_-1^+1 dx
kind=1;
nipe=floor((2*(orderH1+1)+1)/2)+1;
if nipe>20
    error(message('increase nipe in pec.m'));
end
kpts=0;
endpts(1)=-1;
endpts(2)=1;
[bb,xie,w]=gaussq1(kind,nipe,0,0,kpts,endpts);

% generate basis functions on reference edges for applying Dirichlet BCs
[ephdx1,ephdy1,ephdz1,ephdx2,ephdy2,ephdz2,ephdfx1,ephdfy1,ephdfz1,ephdfx2,ephdfy2,ephdfz2,...
    phh1f1,phh1f2,gphfx1,gphfy1,gphfz1,gphfx2,gphfy2,gphfz2,phH1e1,phH1e2,phH1f1,phH1f2]=direvaluate(nipe,xie,order,orderH1,...
    esizet,esizeH1,nipf,intfxi,intfet,gorder);


%--------------------------------------------------------------------------
% Determine vertex based Dirichlet contribution values
%--------------------------------------------------------------------------

fun=probdata.es.dirfun;
arg=probdata.es.dirfunarg;

for i=1:nNodes
        for k=1:3                        % 3 dimensions
            if unkvertexMech(i,k)<0      % This is a Dirichlet DOF
                value=fun(coord(i,1),coord(i,2),coord(i,3),1,arg,probstatic,probFlag);
                dirval(abs(unkvertexMech(i,k)))=value(k);
            end
        end
    
end

%--------------------------------------------------------------------------
% Find edge coefficents
%--------------------------------------------------------------------------
for i=1:nelem
    if orderH1>=2
    
    % transfer coordinates
    xy = coord(intma(i,1:4),1:3);
    
    flag=0;
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
    
    gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
    
    mycoord = zeros(gesizet,3);
    
    %transfer coefficents to locations vertices
    mycoord(1:4,1:3) = xy(1:4,1:3);
    
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
    
    
            
    
    for j=1:6
        if bcedge(glob(i,j))==2
            if j==1
                if eltype(i)==1
                    xi1=1;
                    et1=0;
                    zt1=0;
                    xi2=0;
                    et2=sqrt(3);
                    zt2=0;
                else
                    xi2=1;
                    et2=0;
                    zt2=0;
                    xi1=0;
                    et1=sqrt(3);
                    zt1=0;
                end
            elseif j==2
                xi1=-1;
                et1=0;
                zt1=0;
                xi2=0;
                et2=sqrt(3);
                zt2=0;
            elseif j==3
                xi1=-1;
                et1=0;
                zt1=0;
                xi2=1;
                et2=0;
                zt2=0;
            elseif j==4
                xi1=-1;
                et1=0;
                zt1=0;
                xi2=0;
                et2=sqrt(3)/3;
                zt2=2*(sqrt(2)/sqrt(3));
            elseif j==5
                xi1=1;
                et1=0;
                zt1=0;
                xi2=0;
                et2=sqrt(3)/3;
                zt2=2*(sqrt(2)/sqrt(3));
            else
                xi1=0;
                et1=sqrt(3);
                zt1=0;
                xi2=0;
                et2=sqrt(3)/3;
                zt2=2.*(sqrt(2)/sqrt(3));
            end
            
            % zero array
            al = zeros(orderH1-1,orderH1-1);
            rl = zeros(orderH1-1,3);
            
            
            % apply boundary condition on this edge
            for p=1:nipe
                
                if j==1
                    xi=(0.5*(1-xie(p))*1)+(0.5*(1+xie(p))*0);
                    eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*sqrt(3));
                    zeta=0;
                elseif j==2
                    xi=(0.5*(1.-xie(p))*(-1))+(0.5*(1+xie(p))*0);
                    eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*sqrt(3));
                    zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*0);
                elseif j==3
                    xi=(0.5*(1.-xie(p))*(-1))+(0.5*(1+xie(p))*1);
                    eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*0);
                    zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*0);
                elseif j==4
                    xi=(0.5*(1.-xie(p))*(-1))+(0.5*(1+xie(p))*0);
                    eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*(sqrt(3.)/3));
                    zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*2*(sqrt(2)/sqrt(3)));
                elseif j==5
                    xi=(0.5*(1.-xie(p))*1)+(0.5*(1+xie(p))*0);
                    eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*(sqrt(3)/3));
                    zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*2*(sqrt(2)/sqrt(3)));
                else
                    xi=(0.5*(1.-xie(p))*0)+(0.5*(1+xie(p))*0);
                    eta=(0.5*(1-xie(p))*sqrt(3))+(0.5*(1+xie(p))*(sqrt(3)/3));
                    zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*2*(sqrt(2)/sqrt(3)));
                end
                
                % obtain mapping
                [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian(xy,xi,eta,zeta,...
                    eltype(i),lec,lfc,flag,gorder,mycoord) ;
                
                % compute line det
                ldet=0.5*det*sqrt(...
                    (((axi(1)^2)+(aeta(1)^2)+(azeta(1)^2))*((xi2-xi1)^2))+...
                    (((axi(2)^2)+(aeta(2)^2)+(azeta(2)^2))*((et2-et1)^2))+...
                    (((axi(3)^2)+(aeta(3)^2)+(azeta(3)^2))*((zt2-zt1)^2)));
                
%                  ldet=det*sqrt(...
%                     (((axi(1)^2)+(aeta(1)^2)+(azeta(1)^2))*((xi2-xi1)^2))+...
%                     (((axi(2)^2)+(aeta(2)^2)+(azeta(2)^2))*((et2-et1)^2))+...
%                     (((axi(3)^2)+(aeta(3)^2)+(azeta(3)^2))*((zt2-zt1)^2)));


              
                
                %-----------------------------------------------------------------
                % Obtain H1 basis evaluated on edges
                %-----------------------------------------------------------------
                
                if eltype(i)==1
                   Nh(1:esizeH1,1)=phH1e1((j-1)*nipe+p,1:esizeH1)';
               
                else
                   Nh(1:esizeH1,1)=phH1e2((j-1)*nipe+p,1:esizeH1)';
                end
               
                %-----------------------------------------------------------------
                % Compute vertex contributions on edge
                %-----------------------------------------------------------------
                
                vlow=zeros(3,1);
                
                for pp=1:4
                    for k=1:3
                        if unkvertexMech((intma(i,pp)),k)<0
                            vlow(k)=vlow(k)+(Nh(pp)*dirval(abs(unkvertexMech((intma(i,pp)),k))));
                        end
                    end
                end
          
                
                % computx x,y,z
                [x,y,z]= getxyzcu(xy,xi,eta,zeta,lec,lfc,flag,gorder,eltype(i));
                
             %----------------------------------------------------------------------
             % Compute Dirichlet values from function
             %------------------------------------------------------------------------
                fun=probdata.es.dirfun;
                arg=probdata.es.dirfunarg;
                index=bcedge(glob(i,j));
                value=fun(x,y,z,index,arg,probstatic,probFlag);
                
                
            
                
                %----------------------------------------------------------------------
                % Define system matrix and rhs
                %-----------------------------------------------------------------------
                
                for pp=1:orderH1-1
                    for ppp=1:orderH1-1
                        al(pp,ppp)=al(pp,ppp)+(w(p)*ldet*Nh(j+4+6*(pp-1))*Nh(j+4+6*(ppp-1)));
                    end
            
                    for k=1:3
                        rl(pp,k)=rl(pp,k)-(w(p)*ldet*Nh(j+4+6*(pp-1))*(vlow(k)-value(k)));
                    end
               
                end

            end   
            %----------------------------------------------------------------------------
            % Solve system to get edge coefficients
            %----------------------------------------------------------------------------
            sl = al\rl;
          
            
            %----------------------------------------------------------------------------
            % Put coefficents in matrix system
            %----------------------------------------------------------------------------
          
           

             for pp=1:orderH1-1
                 if abs(unkedgesMech(glob(i,j),pp))~=0
                     %if unkedgesMech(glob(i,j),pp)<0
                    dirval(abs(unkedgesMech(glob(i,j),pp,:)))=sl(pp,:);
                  
                 end
             end
            
        
    end
    end
    end
end


display('corrected edge basis');

% now correct the edge-face basis functions
if orderH1>=3
    
       % normals on refernce element
    nmf(1,1)=-0.5;
    nmf(1,2)=-(sqrt(3)/6);
    nmf(1,3)=-(sqrt(6)/12);
    
    nmf(2,1)=0.5;
    nmf(2,2)=-(sqrt(3)/6);
    nmf(2,3)=-(sqrt(6)/12);
    
    nmf(3,1)=0;
    nmf(3,2)=(sqrt(3)/3);
    nmf(3,3)=-(sqrt(6)/12);
    
    nmf(4,1)=0;
    nmf(4,2)=0;
    nmf(4,3)=(sqrt(3)/(2*sqrt(2)));
    
    % vertices on reference element
    v(1,1)=-1;
    v(1,2)=0;
    v(1,3)=0;
    
    v(2,1)=1;
    v(2,2)=0;
    v(2,3)=0;
    
    v(3,1)=0;
    v(3,2)=sqrt(3);
    v(3,3)=0;
    
    v(4,1)=0;
    v(4,2)=sqrt(3)/3;
    v(4,3)=2*(sqrt(2)/sqrt(3));
    
    for i=1:nelem
        
        % transfer coordinates
        xy = coord(intma(i,1:4),1:3);
        
        % blending coefficents
        flag=0;
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
        
        
        gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
        
        mycoord = zeros(gesizet,3);
        
        %transfer coefficents to locations vertices
        mycoord(1:4,1:3) = xy(1:4,1:3);
        
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
        
        
        if flag==0
            gesizet=4;
        end
        
        if eltype(i)==1
            gphfx=gphfx1;
            gphfy=gphfy1;
            gphfz=gphfz1;
            phh1f=phh1f1;
        else
            gphfx=gphfx2;
            gphfy=gphfy2;
            gphfz=gphfz2;
            phh1f=phh1f2;
        end
        
        for j=1:4
            if cond(i,j)==2
                % perform face correction
                % nf contains local face number
                
                % zero matrix
                alf = zeros(((orderH1-1)^2-(orderH1-1))/2,((orderH1-1)^2-(orderH1-1))/2);
                rlf = zeros(((orderH1-1)^2-(orderH1-1))/2,3);
               
                
                if j==1
                    % set up integration point locations
                    pt(1,1:3)=v(2,1:3);
                    pt(2,1:3)=v(3,1:3);
                    pt(3,1:3)=v(4,1:3);
                    
                elseif j==2
                    pt(1,1:3)=v(3,1:3);
                    pt(2,1:3)=v(1,1:3);
                    pt(3,1:3)=v(4,1:3);
                    
                elseif j==3
                    pt(1,1:3)=v(1,1:3);
                    pt(2,1:3)=v(2,1:3);
                    pt(3,1:3)=v(4,1:3);
                    
                else
                    pt(1,1:3)=v(1,1:3);
                    pt(2,1:3)=v(3,1:3);
                    pt(3,1:3)=v(2,1:3);
                end
                
                if j==1
                    v1(1:3)=xy(3,1:3)-xy(2,1:3);
                    v2(1:3)=xy(4,1:3)-xy(2,1:3);
                    
                elseif j==2
                    v1(1:3)=xy(3,1:3)-xy(1,1:3);
                    v2(1:3)=xy(4,1:3)-xy(1,1:3);
                    
                elseif j==3
                    v1(1:3)=xy(1,1:3)-xy(2,1:3);
                    v2(1:3)=xy(4,1:3)-xy(2,1:3);
                    
                else
                    v2(1:3)=xy(2,1:3)-xy(3,1:3);
                    v1(1:3)=xy(1,1:3)-xy(3,1:3);
                end
                
                for pp=1:nipf
                    % write(6,*)intfxi(pp),intfet(pp),intfw(pp)
                    % compute integration point locations
                    l(1)=1-intfxi(pp)-intfet(pp);
                    l(2)=intfxi(pp);
                    l(3)=intfet(pp);
                    
                    xi = l(1:3)*pt(1:3,1);
                    eta = l(1:3)*pt(1:3,2);
                    zeta = l(1:3)*pt(1:3,3);
                    
                    gph(1:gesizet,1)=gphfx((j-1)*nipf+pp,1:gesizet)';
                    gph(1:gesizet,2)=gphfy((j-1)*nipf+pp,1:gesizet)';
                    gph(1:gesizet,3)=gphfz((j-1)*nipf+pp,1:gesizet)';
                    ph1(1:gesizet,1)=phh1f((j-1)*nipf+pp,1:gesizet)';
                    
                    
                    % evaluate covairant mapping
                    %                     [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian(xy,xi,eta,zeta,...
                    %                         eltype(i),lec,lfc,flag,gorder,mycoord) ;
                    [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord) ;
                    
                    
                    %------------------------------------------------------------------------------------------
                    % Define basis on the face
                    %------------------------------------------------------------------------------------------
                    
                    if eltype(i)==1
                      Nhf=phH1f1((j-1)*nipf+pp,1:esizeH1);
                    else
                      Nhf=phH1f2((j-1)*nipf+pp,1:esizeH1);
                    end
            
                    % computx x,y,z
                    [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
                    
                    % compute curvilinear normal
                    % We use the normal to calculate the area
                    nm(1)=(nmf(j,1)*axi(1))+(nmf(j,2)*aeta(1))+(nmf(j,3)*azeta(1));
                    nm(2)=(nmf(j,1)*axi(2))+(nmf(j,2)*aeta(2))+(nmf(j,3)*azeta(2));
                    nm(3)=(nmf(j,1)*axi(3))+(nmf(j,2)*aeta(3))+(nmf(j,3)*azeta(3));
                    
                    nd=sqrt((nm(1)^2)+(nm(2)^2)+(nm(3)^2));
                    nm(1)=nm(1)/nd;
                    nm(2)=nm(2)/nd;
                    nm(3)=nm(3)/nd;
                    
              %----------------------------------------------------------------------------------------------------
              % Compute Dirichlet values from function
              %----------------------------------------------------------------------------------------------------
              
                    fun=probdata.es.dirfun;
                    arg=probdata.es.dirfunarg;
                    index=cond(i,j);
                    value=fun(x,y,z,index,arg,probstatic,probFlag);
                    
                    
                  
             %-----------------------------------------------------------------------------------------------------
             % Define are for integration purposes
             %-----------------------------------------------------------------------------------------------------
                    area=abs(0.5*scvectpd(nm,v1,v2));
                    
             %-----------------------------------------------------------------------------------------------------
             % Compute the system matrix and rhs to compute the face
             % coefficients
             %-----------------------------------------------------------------------------------------------------
                    
                    ik=0;
                    % first row
                    rk=0;
                    for riif=0:orderH1-3
                        for rjjf=0:orderH1-3
                            if riif+rjjf<=orderH1-3
                                ik=ik+1;
                                rk=rk+1;
                                iif=0;
                                row=4+6*(orderH1-1)+((j-1)*((orderH1-1)^2-(orderH1-1))/2)+rk;
                                
                                
                                % type 1 * 1
                                ck=0;
                                for ciif=0:orderH1-3
                                    for cjjf=0:orderH1-3
                                        if ciif+cjjf<=orderH1-3
                                            iif=iif+1;
                                            ck=ck+1;
                                            col=4+6*(orderH1-1)+((j-1)*((orderH1-1)^2-(orderH1-1))/2)+ck;
                                            
                                            alf(ik,iif)=alf(ik,iif)+Nhf(row)*Nhf(col)*intfw(pp)*area;
                                        end
                                    end
                                end
                            end
                        end
                    end
                  
                  %------------------------------------------------------------------------------------------------
                  % Define vertex contribuitions on faces
                     vlowFace=zeros(3,1);
                
                for ppp=1:4
                    for k=1:3
                        if unkvertexMech((intma(i,ppp)),k)<0
                            vlowFace(k)=vlowFace(k)+(Nhf(ppp)*dirval(abs(unkvertexMech((intma(i,ppp)),k))));
                        end
                    end
                end              
       
                    
                
                  %------------------------------------------------------------------------------------------------
                  % Define edge contributions on faces
                  %------------------------------------------------------------------------------------------------
                  
                     ate=zeros(3,1);
                        for q=1:6
                         
                            for qq=1:orderH1-1
                                for k=1:3
                                if unkedgesMech(glob(i,q),qq)<0
                                    ate(k)=ate(k)+Nhf(4+q+6*(qq-1))*dirval(abs(unkedgesMech(glob(i,q),qq,k)));
                                end
                                end
                            end
                        end
                        
                
                        
                  %------------------------------------------------------------------------------------------------
                  % Define rhs
                  %------------------------------------------------------------------------------------------------
                  
                        ik=0;
                        ck=0;
                        for ciif=0:orderH1-3
                            for cjjf=0:orderH1-3
                                if ciif+cjjf<=orderH1-3
                                    ck=ck+1;
                                    ik=ik+1;
                                    row=4+6*(orderH1-1)+((j-1)*((orderH1-1)^2-(orderH1-1))/2)+ck;
                                    for k=1:3
                                    rlf(ik,k)=rlf(ik,k)-(Nhf(row)*intfw(pp)*area)*(vlowFace(k)+ate(k)-value(k));
                                    end
                                end
                            end
                        end
                  
                    
                end
                
           %----------------------------------------------------------------------------------------------------
           % Solve system to get face coefficients
           %----------------------------------------------------------------------------------------------------
           
                slf = alf\rlf;
             
                
      
           %----------------------------------------------------------------------------------------------------
           % Put coefficents in matrix system
           %----------------------------------------------------------------------------------------------------
           
                ik=0;
                ck=0;
                for ciif=0:orderH1-3
                    for cjjf=0:orderH1-3
                        if ciif+cjjf<=orderH1-3
                            ck=ck+1;
                            ik=ik+1;
                            row=unkfacesMech(globfa(i,j),ck,:);
                         
                            
                                if abs(unkfacesMech(globfa(i,j),ck))~=0
                               %if unkfacesMech(globfa(i,j),ck)<0
                                    dirval(abs(row))=slf(ik,:);
                                    %dirval(abs(row))=0;
                                end
                            
                        end
                    end
                end
                
             
                
              
            end
        end
    end
end

display('completed computation of Dirichlet DOFs')

% change magnitude of diagonal term


DirichletFaceBasis.ephdfx1=ephdfx1;
DirichletFaceBasis.ephdfy1=ephdfy1;
DirichletFaceBasis.ephdfz1=ephdfz1;
DirichletFaceBasis.ephdfx2=ephdfx2;
DirichletFaceBasis.ephdfy2=ephdfy2;
DirichletFaceBasis.ephdfz2=ephdfz2;
DirichletFaceBasis.phh1f1=phh1f1;
DirichletFaceBasis.phh1f2=phh1f2;
DirichletFaceBasis.gphfx1=gphfx1;
DirichletFaceBasis.gphfy1=gphfy1;
DirichletFaceBasis.gphfz1=gphfz1;
DirichletFaceBasis.gphfx2=gphfx2;
DirichletFaceBasis.gphfy2=gphfy2;
DirichletFaceBasis.gphfz2=gphfz2;

% do i=1,npec
% write(6,*)(known(i,j),j=1,nrhs)
% enddo
% pause