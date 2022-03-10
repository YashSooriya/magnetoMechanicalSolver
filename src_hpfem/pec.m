function [known,DirichletFaceBasis]= DirichletEM(mesh,unknown,Quadrature,gorder,order,orderH1,esizet,esizeH1,nipf,nrhs,probdata,probstatic)


% Extract relevant data form mesh structure
nelem=mesh.Nelements;
nside=mesh.edge.nside;
coord=mesh.Coordinates;
glob=mesh.edge.glob;
globfa=mesh.face.globfa;
intma=mesh.intma;
edgecof=mesh.edgecof;
facecof=mesh.facecof;
eltype=mesh.eltype;
bcedge=mesh.edge.bcedge;
cond=mesh.face.cond;

% Extract relevant data from unknown structure
unksid=unknown.unksid;
npec=unknown.npec;
unkfatp1=unknown.unkfatp1;
unkfatp2=unknown.unkfatp2;
unkfatp3=unknown.unkfatp3;
unkz=unknown.unkz;

% Extract relevant data from Quadrature structure
intfxi=Quadrature.intfxi;
intfet=Quadrature.intfet;
intfw=Quadrature.intfw;

probFlag=probdata.probFlag;

lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);
known = zeros(npec,nrhs);

% find integration points along edge
% int_-1^+1 dx
kind=1;
nipe=floor((2*(order+1)+1)/2)+1;
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

% find edge coefficents
for i=1:nelem
    
    if eltype(i)==1
        nme(1,1)=-(1/2);
        nme(1,2)=(sqrt(3))/2;
        nme(1,3)=0;
    else
        nme(1,1)=1/2;
        nme(1,2)=-(sqrt(3))/2;
        nme(1,3)=0;
    end
    nme(2,1)=0.5;
    nme(2,2)=sqrt(3)/2;
    nme(2,3)=0;
    
    nme(3,1)=1;
    nme(3,2)=0;
    nme(3,3)=0;
    
    nme(4,1)=1/2;
    nme(4,2)=sqrt(3)/6;
    nme(4,3)=sqrt(2)/sqrt(3);
    
    nme(5,1)=-1/2;
    nme(5,2)=sqrt(3)/6;
    nme(5,3)=sqrt(2)/sqrt(3);
    
    nme(6,1)=0;
    nme(6,2)=-(2*sqrt(3))/6;
    nme(6,3)=(sqrt(2)/sqrt(3));
    
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
            al = zeros(order+1,order+1);
            rl = zeros(order+1,nrhs);
            
            
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
                    (((asxi(1)^2)+(aseta(1)^2)+(aszeta(1)^2))*((xi2-xi1)^2))+...
                    (((asxi(2)^2)+(aseta(2)^2)+(aszeta(2)^2))*((et2-et1)^2))+...
                    (((asxi(3)^2)+(aseta(3)^2)+(aszeta(3)^2))*((zt2-zt1)^2)));
                % ldet=1.
                
                % obtain basis
                %                ph= basis(xi,eta,zeta,axi,aeta,azeta,order,eltype(i),esizet);
                if eltype(i)==1
                    ph(1:esizet,1:3)=(ephdx1((j-1)*nipe+p,1:esizet)'*axi(1:3))+...
                        (ephdy1((j-1)*nipe+p,1:esizet)'*aeta(1:3))+...
                        (ephdz1((j-1)*nipe+p,1:esizet)'*azeta(1:3));
                else
                    ph(1:esizet,1:3)=(ephdx2((j-1)*nipe+p,1:esizet)'*axi(1:3))+...
                        (ephdy2((j-1)*nipe+p,1:esizet)'*aeta(1:3))+...
                        (ephdz2((j-1)*nipe+p,1:esizet)'*azeta(1:3));
                end
                
                % obtain curvilinear tangent
                nm(1)=(nme(j,1)*asxi(1))+(nme(j,2)*aseta(1))+(nme(j,3)*aszeta(1));
                nm(2)=(nme(j,1)*asxi(2))+(nme(j,2)*aseta(2))+(nme(j,3)*aszeta(2));
                nm(3)=(nme(j,1)*asxi(3))+(nme(j,2)*aseta(3))+(nme(j,3)*aszeta(3));
                
                nd=sqrt((nm(1)^2)+(nm(2)^2)+(nm(3)^2));
                nm(1)=nm(1)/nd;
                nm(2)=nm(2)/nd;
                nm(3)=nm(3)/nd;
                
                % computx x,y,z
                [x,y,z]= getxyzcu(xy,xi,eta,zeta,lec,lfc,flag,gorder,eltype(i));
                
                % compute electric field
                %  use the problem file
                fun=probdata.es.dirfun;
                arg=probdata.es.dirfunarg;
                index=bcedge(glob(i,j));
                e=fun(x,y,z,index,arg,probstatic,probFlag);
                
                
                % tangential exact over an edge
                %te(1,1:nrhs)=nm(1)*e(1)*ones(1,nrhs);
                %te(2,1:nrhs)=nm(2)*e(2)*ones(1,nrhs);
                %te(3,1:nrhs)=nm(3)*e(3)*ones(1,nrhs);
                te(1,1:nrhs)=nm(1)*e(1,1:nrhs);
                te(2,1:nrhs)=nm(2)*e(2,1:nrhs);
                te(3,1:nrhs)=nm(3)*e(3,1:nrhs);
                
                % set up arrays
                for q=0:order
                    for qq=0:order
                        al(q+1,qq+1)=al(q+1,qq+1)+...
                            (((nm(1)*ph(j+(6*q),1))+(nm(2)*ph(j+(6*q),2))+...
                            (nm(3)*ph(j+(6*q),3)))*...
                            ((nm(1)*ph(j+(6*qq),1))+(nm(2)*ph(j+(6*qq),2))+...
                            (nm(3)*ph(j+(6*qq),3)))*w(p)*ldet);
                    end
                end
                % rhs vector
                for iang=1:nrhs
                    for q=0:order
                        rl(q+1,iang)=rl(q+1,iang)+(((nm(1)*ph(j+(6*q),1))+...
                            (nm(2)*ph(j+(6*q),2))+...
                            (nm(3)*ph(j+(6*q),3)))*(te(1,iang)+te(2,iang)+...
                            te(3,iang))*w(p)*ldet);
                    end
                end
                
            end
            sl = al\rl;
            
            % put coefficents in matrix system:-
            
            for q=0:order
                if q==0
                    row=unkz(glob(i,j));
                else
                    row=unksid(glob(i,j),q);
                end
                if row>0
                    error(message('error finding known values'));
                end
                
                % it is possible that all orders may not be included
                for iang=1:nrhs
                    if row~=0
                        known(abs(row),iang)=sl(q+1,iang);
                    end
                end
            end
            
        end
    end
end

display('corrected edge basis');

% now correct the edge-face basis functions
if order>=2
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
                alf = zeros(3*(order-1)+((order-1)*(order-2)),3*(order-1)+((order-1)*(order-2)));
                rlf = zeros(3*(order-1)+((order-1)*(order-2)),nrhs);
                
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
                    
                    
                    
                    % evaluate basis
                    %                    ph= basis(xi,eta,zeta,axi,aeta,azeta,order,eltype(i),esizet);
                    if eltype(i)==1
                        ph(1:esizet,1:3)=(ephdfx1((j-1)*nipf+pp,1:esizet)'*axi(1:3))+...
                            (ephdfy1((j-1)*nipf+pp,1:esizet)'*aeta(1:3))+...
                            (ephdfz1((j-1)*nipf+pp,1:esizet)'*azeta(1:3));
                    else
                        ph(1:esizet,1:3)=(ephdfx2((j-1)*nipf+pp,1:esizet)'*axi(1:3))+...
                            (ephdfy2((j-1)*nipf+pp,1:esizet)'*aeta(1:3))+...
                            (ephdfz2((j-1)*nipf+pp,1:esizet)'*azeta(1:3));
                    end
                    %
                    % computx x,y,z
                    %[x,y,z]= getxyzcu(xy,xi,eta,zeta,lec,lfc,flag,gorder,eltype(i));
                    [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
                    
                    % compute curvilinear normal
                    nm(1)=(nmf(j,1)*axi(1))+(nmf(j,2)*aeta(1))+(nmf(j,3)*azeta(1));
                    nm(2)=(nmf(j,1)*axi(2))+(nmf(j,2)*aeta(2))+(nmf(j,3)*azeta(2));
                    nm(3)=(nmf(j,1)*axi(3))+(nmf(j,2)*aeta(3))+(nmf(j,3)*azeta(3));
                    
                    nd=sqrt((nm(1)^2)+(nm(2)^2)+(nm(3)^2));
                    nm(1)=nm(1)/nd;
                    nm(2)=nm(2)/nd;
                    nm(3)=nm(3)/nd;
                    
                    % compute electric field
                    %  use the problem file
                    fun=probdata.es.dirfun;
                    arg=probdata.es.dirfunarg;
                    index=cond(i,j);
                    e=fun(x,y,z,index,arg,probstatic,probFlag);
                    
                    
                    % tangential exact over face
                    %te(1,1:nrhs)=((nm(2)*e(3))-(nm(3)*e(2)))*ones(1,nrhs);
                    %te(2,1:nrhs)=((nm(3)*e(1))-(nm(1)*e(3)))*ones(1,nrhs);
                    %te(3,1:nrhs)=((nm(1)*e(2))-(nm(2)*e(1)))*ones(1,nrhs);
                    te(1,1:nrhs)=((nm(2)*e(3,1:nrhs))-(nm(3)*e(2,1:nrhs)));
                    te(2,1:nrhs)=((nm(3)*e(1,1:nrhs))-(nm(1)*e(3,1:nrhs)));
                    te(3,1:nrhs)=((nm(1)*e(2,1:nrhs))-(nm(2)*e(1,1:nrhs)));
                    
                    area=abs(0.5*scvectpd(nm,v1,v2));
                    
                    % call surdet(eltype(i),xi,eta,zeta,j,xy,lec,lfc,gorder,area)
                    % area=1.
                    
                    ik=0;
                    % first row
                    rk=0;
                    for riif=0:order-2
                        for rjjf=0:order-2
                            if riif+rjjf<=order-2
                                ik=ik+1;
                                rk=rk+1;
                                iif=0;
                                row=6*(order+1)+((j-1)*(order*order-order)/2)+rk;
                                
                                % type 1 * 1
                                ck=0;
                                for ciif=0:order-2
                                    for cjjf=0:order-2
                                        if ciif+cjjf<=order-2
                                            iif=iif+1;
                                            ck=ck+1;
                                            col=6*(order+1)+((j-1)*...
                                                (order*order-order)/2)+ck;
                                            
                                            alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                                ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                                (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                                ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                                (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                                ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                        end
                                    end
                                end
                                
                                % type 1* 2
                                ck=0;
                                for ciif=0:order-2
                                    for cjjf=0:order-2
                                        if ciif+cjjf<=order-2
                                            iif=iif+1;
                                            ck=ck+1;
                                            col=6*(order+1)+(4*(order*order-order)/2)+((j-1)*(order*order-order)/2)+ck;
                                            
                                            alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                                ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                                (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                                ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                                (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                                ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                        end
                                    end
                                end
                                
                                % type 1* 3
                                for cjjf=0:order-2
                                    iif=iif+1;
                                    col=6*(order+1)+4*((order*order-order)/2)+(4*(order*order-order)/2)+(j-1)*(order-1)+cjjf+1;
                                    
                                    alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                        ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                        (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                        ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                        (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                        ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                end
                                
                            end
                        end
                    end
                    
                    % second row
                    rk=0;
                    for riif=0:order-2
                        for rjjf=0:order-2
                            if riif+rjjf<=order-2
                                ik=ik+1;
                                rk=rk+1;
                                iif=0;
                                row=6*(order+1)+(4*(order*order-order)/2)+((j-1)*(order*order-order)/2)+rk;
                                
                                % type 2 *1
                                ck=0;
                                for ciif=0:order-2
                                    for cjjf=0:order-2
                                        if ciif+cjjf<=order-2
                                            iif=iif+1;
                                            ck=ck+1;
                                            col=6*(order+1)+((j-1)*(order*order-order)/2)+ck;
                                            
                                            
                                            alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                                ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                                (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                                ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                                (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                                ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                        end
                                    end
                                end
                                
                                % type 2  *2
                                ck=0;
                                for ciif=0:order-2
                                    for cjjf=0:order-2
                                        if ciif+cjjf<=order-2
                                            iif=iif+1;
                                            ck=ck+1;
                                            col=6*(order+1)+(4*(order*order-order)/2)+((j-1)*(order*order-order)/2)+ck;
                                            
                                            alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                                ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                                (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                                ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                                (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                                ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                        end
                                    end
                                end
                                
                                % type 2* 3
                                for cjjf=0:order-2
                                    iif=iif+1;
                                    
                                    col=6*(order+1)+4*((order*order-order)/2)+(4*(order*order-order)/2)+(j-1)*(order-1)+cjjf+1;
                                    
                                    alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                        ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                        (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                        ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                        (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                        ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                end
                                
                            end
                        end
                    end
                    
                    % third row
                    rk=0;
                    for rjjf=0:order-2
                        ik=ik+1;
                        iif=0;
                        row=6*(order+1)+4*((order*order-order)/2)+(4*(order*order-order)/2)+(j-1)*(order-1)+rjjf+1;
                        
                        % type 3 *1
                        ck=0;
                        for ciif=0:order-2
                            for cjjf=0:order-2
                                if ciif+cjjf<=order-2
                                    iif=iif+1;
                                    ck=ck+1;
                                    
                                    col=6*(order+1)+((j-1)*(order*order-order)/2)+ck;
                                    
                                    alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                        ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                        (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                        ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                        (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                        ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                    
                                end
                            end
                        end
                        
                        % type 3  *2
                        ck=0;
                        for ciif=0:order-2
                            for cjjf=0:order-2
                                if ciif+cjjf<=order-2
                                    iif=iif+1;
                                    ck=ck+1;
                                    col=6*(order+1)+(4*(order*order-order)/2)+((j-1)*(order*order-order)/2)+ck;
                                    
                                    
                                    alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                        ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                        (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                        ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                        (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                        ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                                    
                                end
                            end
                        end
                        
                        % type 3* 3
                        for cjjf=0:order-2
                            iif=iif+1;
                            col=6*(order+1)+4*((order*order-order)/2)+(4*(order*order-order)/2)+(j-1)*(order-1)+cjjf+1;
                            
                            alf(ik,iif)=alf(ik,iif)+((((nm(2)*ph(row,3))-(nm(3)*ph(row,2)))*...
                                ((nm(2)*ph(col,3))-(nm(3)*ph(col,2))))+ ...
                                (((nm(3)*ph(row,1))-(nm(1)*ph(row,3)))*...
                                ((nm(3)*ph(col,1))-(nm(1)*ph(col,3))))+ ...
                                (((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                ((nm(1)*ph(col,2))-(nm(2)*ph(col,1)))))*intfw(pp)*area;
                        end
                        
                    end
                    
                    for iang=1:nrhs
                        % approximate tangetial given by edge functions
                        ate(1)=complex(0,0);
                        ate(2)=complex(0,0);
                        ate(3)=complex(0,0);
                        for q=1:6
                            if unkz(glob(i,q))<0
                                ate(1)=ate(1)+(((nm(2)*ph(q,3))-(nm(3)*ph(q,2)))*known(abs(unkz(glob(i,q))),iang));
                                ate(2)=ate(2)+(((nm(3)*ph(q,1))-(nm(1)*ph(q,3)))*known(abs(unkz(glob(i,q))),iang));
                                ate(3)=ate(3)+(((nm(1)*ph(q,2))-(nm(2)*ph(q,1)))*known(abs(unkz(glob(i,q))),iang));
                            end
                            for qq=1:order
                                if unksid(glob(i,q),qq)<0
                                    ate(1)=ate(1)+(((nm(2)*ph(q+(6*qq),3))-(nm(3)*ph(q+(6*qq),2)))*known(abs(unksid(glob(i,q),qq)),iang));
                                    ate(2)=ate(2)+(((nm(3)*ph(q+(6*qq),1))-(nm(1)*ph(q+(6*qq),3)))*known(abs(unksid(glob(i,q),qq)),iang));
                                    ate(3)=ate(3)+(((nm(1)*ph(q+(6*qq),2))-(nm(2)*ph(q+(6*qq),1)))*known(abs(unksid(glob(i,q),qq)),iang));
                                end
                            end
                        end
                        
                        ik=0;
                        ck=0;
                        for ciif=0:order-2
                            for cjjf=0:order-2
                                if ciif+cjjf<=order-2
                                    ck=ck+1;
                                    ik=ik+1;
                                    row=6*(order+1)+((j-1)*(order*order-order)/2)+ck;
                                    
                                    rlf(ik,iang)=rlf(ik,iang)+(((((nm(2)*ph(row ,3))-(nm(3)*ph(row ,2)))*...
                                        (te(1,iang)-ate(1)))+(((nm(3)*ph(row ,1))-(nm(1)*ph(row ,3)))*...
                                        (te(2,iang)-ate(2)))+(((nm(1)*ph(row ,2))-(nm(2)*ph(row ,1)))*...
                                        (te(3,iang)-ate(3))))*intfw(pp)*area);
                                end
                            end
                        end
                        
                        ck=0;
                        for ciif=0:order-2
                            for cjjf=0:order-2
                                if ciif+cjjf<=order-2
                                    ck=ck+1;
                                    ik=ik+1;
                                    row=6*(order+1)+(4*(order*order-order)/2)+((j-1)*(order*order-order)/2)+ck;
                                    
                                    rlf(ik,iang)=rlf(ik,iang)+(((((nm(2)*ph(row ,3))-(nm(3)*ph(row ,2)))*...
                                        (te(1,iang)-ate(1)))+(((nm(3)*ph(row ,1))-(nm(1)*ph(row ,3)))*...
                                        (te(2,iang)-ate(2)))+(((nm(1)*ph(row,2))-(nm(2)*ph(row,1)))*...
                                        (te(3,iang)-ate(3))))*intfw(pp)*area);
                                end
                            end
                        end
                        
                        for cjjf=0:order-2
                            ik=ik+1;
                            row=6*(order+1)+4*((order*order-order)/2)+(4*(order*order-order)/2)+(j-1)*(order-1)+cjjf+1;
                            
                            rlf(ik,iang)=rlf(ik,iang)+(((((nm(2)*ph(row ,3))-(nm(3)*ph(row ,2)))*...
                                (te(1,iang)-ate(1)))+(((nm(3)*ph(row ,1))-(nm(1)*ph(row ,3)))*...
                                (te(2,iang)-ate(2)))+(((nm(1)*ph(row,2))-(nm(2)*ph(row ,1)))*...
                                (te(3,iang)-ate(3))))*intfw(pp)*area);
                            
                        end
                        
                        
                    end
                    
                end
                
                %solve system:-
                slf = alf\rlf;
                
                % write(6,*)i,j
                % write(6,*)(slf(if,1),if=1,(order+1)*(order-1))
                % pause
                
                % put coefficents in matrix system:-
                ik=0;
                ck=0;
                for ciif=0:order-2
                    for cjjf=0:order-2
                        if ciif+cjjf<=order-2
                            ck=ck+1;
                            ik=ik+1;
                            row=unkfatp1(globfa(i,j),ck);
                            if row>0 || abs(row)>npec
                                error(message('error finding known values1'));
                            end
                            % it is possible that all rows might not be included
                            for iang=1:nrhs
                                if row~=0
                                    known(abs(row),iang)=slf(ik,iang);
                                end
                            end
                        end
                    end
                end
                
                ck=0;
                for ciif=0:order-2
                    for cjjf=0:order-2
                        if ciif+cjjf<=order-2
                            ck=ck+1;
                            ik=ik+1;
                            row=unkfatp2(globfa(i,j),ck);
                            if row>0 || abs(row)>npec
                                error(message('error finding known values 2'));
                            end
                            % it is possible that all rows might not be included
                            for iang=1:nrhs
                                % write(6,*)abs(row),iang,known(abs(row),iang),slf(if,iang)
                                if(row~=0)
                                    known(abs(row),iang)=slf(ik,iang);
                                end
                            end
                        end
                    end
                end
                
                ck=0;
                for cjjf=0:order-2
                    ck=ck+1;
                    ik=ik+1;
                    row=unkfatp3(globfa(i,j),ck);
                    if row>0 || abs(row)>npec
                        error(message('error finding known values 3'));
                    end
                    % it is possible that all rows might not be included
                    for iang=1:nrhs
                        if row~=0
                            known(abs(row),iang)=slf(ik,iang);
                        end
                    end
                end
            end
        end
    end
end

display('completed computation of Dirichlet DOFs')

% apply far field condition n x E  = 0
% change magnitude of diagonal term

for i=1:nside
    if bcedge(i)==1
        for j=0:order
            if j==0
                row=unkz(i);
            else
                row=unksid(i,j);
            end
            
            for iang=1:nrhs
                known(abs(row),iang)=complex(0,0);
            end
        end
    end
end

% change magnitude of diagonal term
if order>=2
    for i=1:nelem
        for j=1:4
            if cond(i,j)==1
                for jf=1:(order*order-order)/2
                    row=unkfatp1(globfa(i,j),jf);
                    if row~=0
                        for iang=1:nrhs
                            known(abs(row),iang)=0;
                        end
                    end
                end
                for jf=1:(order*order-order)/2
                    row=unkfatp2(globfa(i,j),jf);
                    if row~=0
                        for iang=1:nrhs
                            known(abs(row),iang)=complex(0,0);
                        end
                    end
                end
                
                for jf=1:order-1
                    row=unkfatp3(globfa(i,j),jf);
                    if row~=0
                        for iang=1:nrhs
                            known(abs(row),iang)=complex(0,0);
                        end
                    end
                end
                
            end
        end
    end
end

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