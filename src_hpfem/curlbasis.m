function ph= curlbasis(xi,eta,zeta,asxi,aseta,aszeta,order,eltype,esize)

ph = zeros(esize,3);
% area coordinates and there derivitives
l(1)=0.5-(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(2)=0.5+(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(3)=((sqrt(3)/3)*eta)-((sqrt(6)/12)*zeta);
l(4)=((sqrt(3)/(2*sqrt(2)))*zeta);

% derivitives of area coordinates
dl(1,1)=-0.5;
dl(1,2)=-sqrt(3)/6;
dl(1,3)=-sqrt(6)/12;

dl(2,1)=0.5;
dl(2,2)=-sqrt(3)/6;
dl(2,3)=-sqrt(6)/12;

dl(3,1)=0;
dl(3,2)=sqrt(3)/3;
dl(3,3)=-sqrt(6)/12;

dl(4,1)=0;
dl(4,2)=0;
dl(4,3)=sqrt(3)/(2*sqrt(2));

% paramertisations
if eltype == 1
    s(1)=l(3)-l(2);
else
    s(1)=l(2)-l(3);
end
s(2)=l(3)-l(1);
s(3)=l(2)-l(1);
s(4)=l(4)-l(1);
s(5)=l(4)-l(2);
s(6)=l(4)-l(3);
	
% derivitives of paramertisations
if eltype==1
    ds(1,1:3)=dl(3,1:3)-dl(2,1:3);
else
    ds(1,1:3)=dl(2,1:3)-dl(3,1:3);
end
ds(2,1:3)=dl(3,1:3)-dl(1,1:3);
ds(3,1:3)=dl(2,1:3)-dl(1,1:3);
ds(4,1:3)=dl(4,1:3)-dl(1,1:3);
ds(5,1:3)=dl(4,1:3)-dl(2,1:3);
ds(6,1:3)=dl(4,1:3)-dl(3,1:3);
      
% p=0 edge basis functions
if eltype==1
    ph(1,1)=((1./sqrt(2))*asxi(1))+((1/sqrt(6))*aseta(1))+((2/sqrt(3))*aszeta(1)); 
    ph(1,2)=((1./sqrt(2))*asxi(2))+((1./sqrt(6))*aseta(2))+((2/sqrt(3))*aszeta(2)); 
    ph(1,3)=((1./sqrt(2))*asxi(3))+((1./sqrt(6))*aseta(3))+((2/sqrt(3))*aszeta(3)); 
else
    ph(1,1)=-((1/sqrt(2))*asxi(1))-((1/sqrt(6))*aseta(1))-((2/sqrt(3))*aszeta(1));
    ph(1,2)=-((1/sqrt(2))*asxi(2))-((1/sqrt(6))*aseta(2))-((2/sqrt(3))*aszeta(2));
    ph(1,3)=-((1/sqrt(2))*asxi(3))-((1/sqrt(6))*aseta(3))-((2/sqrt(3))*aszeta(3)); 
end

ph(2,1)=((1/sqrt(2))*asxi(1))+((-1/sqrt(6))*aseta(1))+((-2/sqrt(3))*aszeta(1));
ph(2,2)=((1/sqrt(2))*asxi(2))+((-1/sqrt(6))*aseta(2))+((-2/sqrt(3))*aszeta(2));
ph(2,3)=((1/sqrt(2))*asxi(3))+((-1/sqrt(6))*aseta(3))+((-2/sqrt(3))*aszeta(3));

ph(3,1)=((-sqrt(2/3))*aseta(1))+((2/sqrt(3))*aszeta(1));
ph(3,2)=((-sqrt(2/3))*aseta(2))+((2/sqrt(3))*aszeta(2));
ph(3,3)=((-sqrt(2/3))*aseta(3))+((2/sqrt(3))*aszeta(3));

ph(4,1)=((-1/sqrt(2))*asxi(1))+((sqrt(3/2))*aseta(1));
ph(4,2)=((-1/sqrt(2))*asxi(2))+((sqrt(3./2.))*aseta(2));
ph(4,3)=((-1./sqrt(2))*asxi(3))+((sqrt(3/2))*aseta(3));

ph(5,1)=((-1/sqrt(2))*asxi(1))+((-sqrt(3/2))*aseta(1));
ph(5,2)=((-1/sqrt(2))*asxi(2))+((-sqrt(3/2))*aseta(2));
ph(5,3)=((-1/sqrt(2))*asxi(3))+((-sqrt(3/2))*aseta(3));

ph(6,1)=((sqrt(2))*asxi(1));
ph(6,2)=((sqrt(2))*asxi(2));
ph(6,3)=((sqrt(2))*asxi(3));
     
% the higher order edge functions have zero curl
% they are gradients of scaler.
% nabla x grad phi =0 for all phi!!!
% higher order edge basis functions
for m=1:6
    for p=1:order        
        ph(m+6*p,1:3)=0;        
    end
end
      
nbas=6*(order+1);
      
if order>=2
    
    % define sfu(m), m=1,2,3,4
    if eltype==1
        sfu(1)=l(2)-l(3);
        sfu(4)=l(1)-l(2);
    else
        sfu(1)=l(3)-l(2);
        sfu(4)=l(1)-l(3);
    end
    sfu(2)=l(1)-l(3);
    sfu(3)=l(1)-l(2);
    
    % define dsfu(m,j), m=1,2,3,4,j=1,2,3
    if eltype==1
        dsfu(1,1:3)=dl(2,1:3)-dl(3,1:3);
        dsfu(4,1:3)=dl(1,1:3)-dl(2,1:3);
    else
        dsfu(1,1:3)=dl(3,1:3)-dl(2,1:3);
        dsfu(4,1:3)=dl(1,1:3)-dl(3,1:3);
    end
    dsfu(2,1:3)=dl(1,1:3)-dl(3,1:3);
    dsfu(3,1:3)=dl(1,1:3)-dl(2,1:3);
    
    
    % define tfu(m), m=1,2,3,4
    if eltype==1
        tfu(1)=l(2)+l(3);
        tfu(4)=l(1)+l(2);
    else
        tfu(1)=l(3)+l(2);
        tfu(4)=l(1)+l(3);
    end
    tfu(2)=l(1)+l(3);
    tfu(3)=l(1)+l(2);
    
    % define dtfu(m,j), m=1,2,3,4,j=1,2,3
    if eltype==1
        dtfu(1,1:3)=dl(2,1:3)+dl(3,1:3);
        dtfu(4,1:3)=dl(1,1:3)+dl(2,1:3);
    else
        dtfu(1,1:3)=dl(3,1:3)+dl(2,1:3);
        dtfu(4,1:3)=dl(1,1:3)+dl(3,1:3);
    end
    dtfu(2,1:3)=dl(1,1:3)+dl(3,1:3);
    dtfu(3,1:3)=dl(1,1:3)+dl(2,1:3);
    
    % define sfv(m), m=1,2,3,4
    if eltype==1
        sfv(1)=l(4)-l(2)-l(3);
        sfv(4)=l(3)-l(2)-l(1);
    else
        sfv(1)=l(4)-l(2)-l(3);
        sfv(4)=l(2)-l(3)-l(1);
    end
    sfv(2)=l(4)-l(3)-l(1);
    sfv(3)=l(4)-l(2)-l(1);
    
    % define dsfv(m,j), m=1,2,3,4,j=1,2,3
    if eltype==1
        dsfv(1,1:3)=dl(4,1:3)-dl(2,1:3)-dl(3,1:3);
        dsfv(4,1:3)=dl(3,1:3)-dl(2,1:3)-dl(1,1:3);
    else
        dsfv(1,1:3)=dl(4,1:3)-dl(2,1:3)-dl(3,1:3);
        dsfv(4,1:3)=dl(2,1:3)-dl(3,1:3)-dl(1,1:3);
    end
    dsfv(2,1:3)=dl(4,1:3)-dl(3,1:3)-dl(1,1:3);
    dsfv(3,1:3)=dl(4,1:3)-dl(2,1:3)-dl(1,1:3);
    
    
    % define tfu(m), m=1,2,3,4
    if eltype==1
        tfv(1)=l(2)+l(3)+l(4);
        tfv(4)=l(1)+l(2)+l(3);
    else
        tfv(1)=l(3)+l(2)+l(4);
        tfv(4)=l(1)+l(3)+l(2);
    end
    tfv(2)=l(1)+l(3)+l(4);
    tfv(3)=l(1)+l(2)+l(4);
    
    % define dtfu(m,j), m=1,2,3,4,j=1,2,3
    if eltype==1
        dtfv(1,1:3)=dl(2,1:3)+dl(3,1:3)+dl(4,1:3);
        dtfv(4,1:3)=dl(1,1:3)+dl(2,1:3)+dl(3,1:3);
    else
        dtfv(1,1:3)=dl(3,1:3)+dl(2,1:3)+dl(4,1:3);
        dtfv(4,1:3)=dl(1,1:3)+dl(3,1:3)+dl(2,1:3);
    end
    dtfv(2,1:3)=dl(1,1:3)+dl(3,1:3)+dl(4,1:3);
    dtfv(3,1:3)=dl(1,1:3)+dl(2,1:3)+dl(4,1:3);
    
    % define lf(m),m=1,2,3,4
    lf(1)=l(4);
    lf(2)=l(4);
    lf(3)=l(4);
    if eltype==1
        lf(4)=l(3);
    else
        lf(4)=l(2);
    end
    
    % define dlf(m,j),m=1,2,3,4,j=1,2,3
    dlf(1,1:3)=dl(4,1:3);
    dlf(2,1:3)=dl(4,1:3);
    dlf(3,1:3)=dl(4,1:3);
    if eltype==1
        dlf(4,1:3)=dl(3,1:3);
    else
        dlf(4,1:3)=dl(2,1:3);
    end
    
    
    % type 1 functions
    % have zero curl (gradients of scalars!)
    nbas=6*(order+1);
    for m=1:4
        for i=0:order-2
            for j=0:order-2
                if i+j<=order-2
                    nbas=nbas+1;
                    ph(nbas,1:3)=0;
                end
            end
        end
    end
    
    %type 2 functions
    for m=1:4
        for i=0:order-2
            for j=0:order-2
                if i+j<=order-2
                    nbas=nbas+1;
                    % define ui
                    if abs(tfu(m))>1e-12
                        ui=iln(i+2,sfu(m)/tfu(m))*(tfu(m)^(i+2));
                    else
                        ui=0;
                    end
                    % define vi
                    if abs(tfv(m))>1e-12
                        vj=lf(m)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j));
                    else
                        vj=0;
                    end
                    % define dui
                    if abs(tfu(m))>1e-12
                        dui(1)=-ln(i,sfu(m)/tfu(m))*(tfu(m)^(i+1))*dtfu(m,1)+...
                            ln(i+1,sfu(m)/tfu(m))*(tfu(m)^(i+1))*dsfu(m,1);
                        dui(2)=-ln(i,sfu(m)/tfu(m))*(tfu(m)^(i+1))*dtfu(m,2)+...
                            ln(i+1,sfu(m)/tfu(m))*(tfu(m)^(i+1))*dsfu(m,2);
                        dui(3)=-ln(i,sfu(m)/tfu(m))*(tfu(m)^(i+1))*dtfu(m,3)+...
                            ln(i+1,sfu(m)/tfu(m))*(tfu(m)^(i+1))*dsfu(m,3);
                    else
                        dui(1)=0;
                        dui(2)=0;
                        dui(3)=0;
                    end
                    % define dvj
                    if abs(tfv(m))>1e-12
                        dvj(1)=dlf(m,1)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                            lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,1)-...
                            sfv(m)*dtfv(m,1))/(tfv(m)^2))*(tfv(m)^j)+...
                            lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,1);
                        
                        dvj(2)=dlf(m,2)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                            lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,2)-...
                            sfv(m)*dtfv(m,2))/(tfv(m)^2))*(tfv(m)^j)+...
                            lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,2);
                        
                        dvj(3)=dlf(m,3)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                            lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,3)-...
                            sfv(m)*dtfv(m,3))/(tfv(m)^2))*(tfv(m)^j)+...
                            lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,3);
                    else
                        dvj(1)=0;
                        dvj(2)=0;
                        dvj(3)=0;
                    end
                    
                    % much simpler calculation!!!
                    dphizy=dui(3)*dvj(2) - ( dvj(3)*dui(2) );
                    dphiyz=dui(2)*dvj(3) - ( dvj(2)*dui(3));
                    
                    dphixz=dui(1)*dvj(3) - ( dvj(1)*dui(3) );
                    dphizx=dui(3)*dvj(1) - ( dvj(3)*dui(1) );
                    
                    dphiyx=dui(2)*dvj(1) - ( dvj(2)*dui(1) );
                    dphixy=dui(1)*dvj(2) - ( dvj(1)*dui(2) );
                    
                    
                    % define derivatives of nablaui
                    % nb dtfu and dsfu are constant vectors
                    
                    ph(nbas,1:3)=(dphizy-dphiyz)*asxi(1:3)+...
                                 (dphixz-dphizx)*aseta(1:3)+...
                                 (dphiyx-dphixy)*aszeta(1:3);
                    
                end
            end
        end
    end
    
    
    % type 3 functions
    % we need to know lft3(1),lft3(2) for each face and their derivatives
    if eltype==1
        lft3(1,1)=l(2);
        lft3(1,2)=l(3);
        lft3(4,1)=l(1);
        lft3(4,2)=l(2);
    else
        lft3(1,1)=l(3);
        lft3(1,2)=l(2);
        lft3(4,1)=l(1);
        lft3(4,2)=l(3);
    end
    lft3(2,1)=l(1);
    lft3(2,2)=l(3);
    lft3(3,1)=l(1);
    lft3(3,2)=l(2);
    % their derivatives are
    if eltype==1
        d1lft3(1,1:3)=dl(2,1:3);
        d2lft3(1,1:3)=dl(3,1:3);
        d1lft3(4,1:3)=dl(1,1:3);
        d2lft3(4,1:3)=dl(2,1:3);
    else
        d1lft3(1,1:3)=dl(3,1:3);
        d2lft3(1,1:3)=dl(2,1:3);
        d1lft3(4,1:3)=dl(1,1:3);
        d2lft3(4,1:3)=dl(3,1:3);
    end
    d1lft3(2,1:3)=dl(1,1:3);
    d2lft3(2,1:3)=dl(3,1:3);
    d1lft3(3,1:3)=dl(1,1:3);
    d2lft3(3,1:3)=dl(2,1:3);
    
    
    for m=1:4
        % simple caluclation added
        if eltype==1
            % type 1
            if m==1
                gf1(1:3)=dl(2,1:3);
                gf2(1:3)=dl(3,1:3);
                
                f1=l(2);
                f2=l(3);
            elseif m==2
                gf1(1:3)=dl(1,1:3);
                gf2(1:3)=dl(3,1:3);
                
                f1=l(1);
                f2=l(3);
            elseif m==3
                gf1(1:3)=dl(1,1:3);
                gf2(1:3)=dl(2,1:3);
                
                f1=l(1);
                f2=l(2);
            else
                gf1(1:3)=dl(1,1:3);
                gf2(1:3)=dl(2,1:3);
                
                f1=l(1);
                f2=l(2);
            end
        else
            % type 2
            if m==1
                gf1(1:3)=dl(3,1:3);
                gf2(1:3)=dl(2,1:3);
                
                f1=l(3);
                f2=l(2);
            elseif m==2
                gf1(1:3)=dl(1,1:3);
                gf2(1:3)=dl(3,1:3);
                
                f1=l(1);
                f2=l(3);
            elseif m==3
                gf1(1:3)=dl(1,1:3);
                gf2(1:3)=dl(2,1:3);
                
                f1=l(1);
                f2=l(2);
            else
                gf1(1:3)=dl(1,1:3);
                gf2(1:3)=dl(3,1:3);
                
                f1=l(1);
                f2=l(3);
            end
        end
        
        
        
        
        for j=0:order-2
            nbas=nbas+1;
            if abs(tfv(m))>1e-12
                % define vi
                vj=lf(m)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j));
                
                % define dvj
                dvj(1)=dlf(m,1)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                    lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,1)-...
                    sfv(m)*dtfv(m,1))/(tfv(m)^2))*(tfv(m)^j)+...
                    lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,1);
                
                dvj(2)=dlf(m,2)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                    lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,2)-...
                    sfv(m)*dtfv(m,2))/(tfv(m)^2))*(tfv(m)^j)+...
                    lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,2);
                
                dvj(3)=dlf(m,3)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                    lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,3)-...
                    sfv(m)*dtfv(m,3))/(tfv(m)^2))*(tfv(m)^j)+...
                    lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,3);
            else
                vj=0;
                dvj(1)=0;
                dvj(2)=0;
                dvj(3)=0;
            end
            
            % compute components of curl
            % much simpler calculation!!!
            dphizy=gf1(3)*(dvj(2)*f2+vj*gf2(2)) - (gf2(3)*(dvj(2)*f1+vj*gf1(2) ));
            dphiyz=gf1(2)*(dvj(3)*f2+vj*gf2(3)) -(gf2(2)*(dvj(3)*f1+vj*gf1(3) ));
            dphixz=gf1(1)*(dvj(3)*f2+vj*gf2(3)) - (gf2(1)*(dvj(3)*f1+vj*gf1(3) ));
            dphizx=gf1(3)*(dvj(1)*f2+vj*gf2(1)) - (gf2(3)*(dvj(1)*f1+vj*gf1(1) ));
            dphiyx=gf1(2)*(dvj(1)*f2+vj*gf2(1)) - (gf2(2)*(dvj(1)*f1+vj*gf1(1) ));
            dphixy=gf1(1)*(dvj(2)*f2+vj*gf2(2)) - (gf2(1)*(dvj(2)*f1+vj*gf1(2) ));
            
            
            
            ph(nbas,1:3)=(dphizy-dphiyz)*asxi(1:3)+...
                         (dphixz-dphizx)*aseta(1:3)+ (dphiyx-dphixy)*aszeta(1:3);
            
        end
    end
end



if order>=3
    si(1)=l(1)-l(2);
    ti(1)=l(1)+l(2);
    si(2)=l(3)-l(1)-l(2);
    ti(2)=l(3)+l(1)+l(2);
    si(3)=l(4)-l(1)-l(2)-l(3);
    
    dsi(1,1)=dl(1,1)-dl(2,1);
    dsi(1,2)=dl(1,2)-dl(2,2);
    dsi(1,3)=dl(1,3)-dl(2,3);
    dti(1,1)=dl(1,1)+dl(2,1);
    dti(1,2)=dl(1,2)+dl(2,2);
    dti(1,3)=dl(1,3)+dl(2,3);
    
    dsi(2,1)=dl(3,1)-dl(1,1)-dl(2,1);
    dsi(2,2)=dl(3,2)-dl(1,2)-dl(2,2);
    dsi(2,3)=dl(3,3)-dl(1,3)-dl(2,3);
    dti(2,1)=dl(3,1)+dl(1,1)+dl(2,1);
    dti(2,2)=dl(3,2)+dl(1,2)+dl(2,2);
    dti(2,3)=dl(3,3)+dl(1,3)+dl(2,3);
    
    dsi(3,1)=dl(4,1)-dl(1,1)-dl(2,1)-dl(3,1);
    dsi(3,2)=dl(4,2)-dl(1,2)-dl(2,2)-dl(3,2);
    dsi(3,3)=dl(4,3)-dl(1,3)-dl(2,3)-dl(3,3);
    
    % gradients
    for i=0:order-3
        for j=0:order-3
            for k=0:order-3
                if i+j+k<=order-3
                    % type 1, zero curl
                    nbas=nbas+1;
                    ph(nbas,1:3)=0;
                end
            end
        end
    end
    
    %non-gradients
    for i=0:order-3
        for j=0:order-3
            for k=0:order-3
                if i+j+k<=order-3
    
                    
                    if abs(ti(1))>1e-12
                        % write(6,*)'type 1',i,j,k
                        ui=iln(i+2,si(1)/ti(1))*(ti(1)^(i+2));
                        dd=((dsi(1,1)*ti(1))-(si(1)*dti(1,1)))/(ti(1)^2);
                        dui(1)=diln(i+2,si(1)/ti(1))*dd*(ti(1)^(i+2))+...
                            (iln(i+2,si(1)/ti(1))*(i+2)*(ti(1)^(i+1))*dti(1,1));
                        dd=((dsi(1,2)*ti(1))-(si(1)*dti(1,2)))/(ti(1)^2);
                        dui(2)=diln(i+2,si(1)/ti(1))*dd*(ti(1)^(i+2))+...
                            (iln(i+2,si(1)/ti(1))*(i+2)*(ti(1)^(i+1))*dti(1,2));
                        dd=((dsi(1,3)*ti(1))-(si(1)*dti(1,3)))/(ti(1)^2);
                        dui(3)=diln(i+2,si(1)/ti(1))*dd*(ti(1)^(i+2))+...
                            (iln(i+2,si(1)/ti(1))*(i+2)*(ti(1)^(i+1))*dti(1,3));
                        % write(6,*)dui(1),dui(2),dui(3)
                        
                        vj=l(3)*ln(j,si(2)/ti(2))*(ti(2)^j);
                        % dd= d/dx(si(2)/ti(2)
                        dd=((dsi(2,1)*ti(2))-(si(2)*dti(2,1)))/(ti(2)^2);
                        dvj(1)=dl(3,1)*ln(j,si(2)/ti(2))*(ti(2)^j)+...
                            l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)+...
                            l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,1);
                        % dd= d/dy(si(2)/ti(2)
                        dd=((dsi(2,2)*ti(2))-(si(2)*dti(2,2)))/(ti(2)^2);
                        dvj(2)=dl(3,2)*ln(j,si(2)/ti(2))*(ti(2)^j)+...
                            l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)+...
                            l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,2);
                        % dd= d/dz(si(2)/ti(2)
                        dd=((dsi(2,3)*ti(2))-(si(2)*dti(2,3)))/(ti(2)^2);
                        dvj(3)=dl(3,3)*ln(j,si(2)/ti(2))*(ti(2)^j)+...
                            l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)+...
                            l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,3);
                    else
                        ui=0;
                        dui(1)=0;
                        dui(2)=0;
                        dui(3)=0;
                        dduixy=0;
                        dduixz=0;
                        dduiyx=0;
                        dduiyz=0;
                        dduizx=0;
                        dduizy=0;
                        
                        vj=0;
                        dvj(1)=0;
                        dvj(2)=0;
                        dvj(3)=0;
                        ddvjxy=0;
                        ddvjxz=0;
                        ddvjyx=0;
                        ddvjyz=0;
                        ddvjzx=0;
                        ddvjzy=0;
                    end
                    
                    wk=l(4)*ln(k,si(3));
                    dwk(1)=dl(4,1)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,1);
                    dwk(2)=dl(4,2)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,2);
                    dwk(3)=dl(4,3)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,3);
                    
                    % much simpler calculation!!!
                    dphizy=dui(3)*(dvj(2)*wk+vj*dwk(2)) - ...
                        (dvj(3)*(dui(2)*wk+ui*dwk(2) )) +...
                        (dwk(3)*(dui(2)*vj+ui*dvj(2) ));
                    
                    dphiyz=dui(2)*(dvj(3)*wk+vj*dwk(3)) - ...
                        (dvj(2)*(dui(3)*wk+ui*dwk(3) )) +...
                        (dwk(2)*(dui(3)*vj+ui*dvj(3) ));
                    
                    dphixz=dui(1)*(dvj(3)*wk+vj*dwk(3)) - ...
                        (dvj(1)*(dui(3)*wk+ui*dwk(3) )) +...
                        (dwk(1)*(dui(3)*vj+ui*dvj(3) ));
                    
                    dphizx=dui(3)*(dvj(1)*wk+vj*dwk(1)) - ...
                        (dvj(3)*(dui(1)*wk+ui*dwk(1) )) +...
                        (dwk(3)*(dui(1)*vj+ui*dvj(1) ));
                    
                    dphiyx=dui(2)*(dvj(1)*wk+vj*dwk(1)) - ...
                        (dvj(2)*(dui(1)*wk+ui*dwk(1) )) +...
                        (dwk(2)*(dui(1)*vj+ui*dvj(1) ));
                    
                    dphixy=dui(1)*(dvj(2)*wk+vj*dwk(2)) - ...
                        (dvj(1)*(dui(2)*wk+ui*dwk(2) )) +...
                        (dwk(1)*(dui(2)*vj+ui*dvj(2) ));
                    
                    % define derivatives of nablaui
                    % nb dtfu and dsfu are constant vectors
                    nbas=nbas+1;
                    ph(nbas,1:3)=(dphizy-dphiyz)*asxi(1:3)+...
                        (dphixz-dphizx)*aseta(1:3)+(dphiyx-dphixy)*aszeta(1:3);
                    
                    % compute components of curl
                    % much simpler calculation!!!
                    dphizy=dui(3)*(dvj(2)*wk+vj*dwk(2)) -...
                        (dvj(3)*(dui(2)*wk+ui*dwk(2) )) -...
                        (dwk(3)*(dui(2)*vj+ui*dvj(2) ));
                    
                    dphiyz=dui(2)*(dvj(3)*wk+vj*dwk(3)) - ...
                        (dvj(2)*(dui(3)*wk+ui*dwk(3) )) -...
                        (dwk(2)*(dui(3)*vj+ui*dvj(3) ));
                    
                    dphixz=dui(1)*(dvj(3)*wk+vj*dwk(3)) - ...
                        (dvj(1)*(dui(3)*wk+ui*dwk(3) )) -...
                        (dwk(1)*(dui(3)*vj+ui*dvj(3) ));
                    
                    dphizx=dui(3)*(dvj(1)*wk+vj*dwk(1)) - ...
                        (dvj(3)*(dui(1)*wk+ui*dwk(1) )) -...
                        (dwk(3)*(dui(1)*vj+ui*dvj(1) ));
                    
                    dphiyx=dui(2)*(dvj(1)*wk+vj*dwk(1)) - ...
                        (dvj(2)*(dui(1)*wk+ui*dwk(1) )) -...
                        (dwk(2)*(dui(1)*vj+ui*dvj(1) ));
                    
                    dphixy=dui(1)*(dvj(2)*wk+vj*dwk(2)) - ...
                        (dvj(1)*(dui(2)*wk+ui*dwk(2) )) -...
                        (dwk(1)*(dui(2)*vj+ui*dvj(2) ));
                    
                    nbas=nbas+1;
                    ph(nbas,1:3)=(dphizy-dphiyz)*asxi(1:3)+...
                        (dphixz-dphizx)*aseta(1:3)+ (dphiyx-dphixy)*aszeta(1:3);
                end
            end
        end
    end
    
    
    for j=0:order-3
        for k=0:order-3
            if j+k<=order-3
                nbas=nbas+1;
                wk=l(4)*ln(k,si(3));
                dwk(1)=dl(4,1)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,1);
                dwk(2)=dl(4,2)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,2);
                dwk(3)=dl(4,3)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,3);
                
                if abs(ti(2))>1e-12
                    vj=l(3)*ln(j,si(2)/ti(2))*(ti(2)^j);
                    dd=((dsi(2,1)*ti(2))-(si(2)*dti(2,1)))/(ti(2)^2);
                    
                    dvj(1)=dl(3,1)*ln(j,si(2)/ti(2))*(ti(2)^j)+...
                        l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)+...
                        l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,1);
                    
                    dd=((dsi(2,2)*ti(2))-(si(2)*dti(2,2)))/(ti(2)^2);
                    
                    dvj(2)=dl(3,2)*ln(j,si(2)/ti(2))*(ti(2)^j)+...
                        l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)+...
                        l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,2);
                    
                    dd=((dsi(2,3)*ti(2))-(si(2)*dti(2,3)))/(ti(2)^2);
                    dvj(3)=dl(3,3)*ln(j,si(2)/ti(2))*(ti(2)^j)+...
                        l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)+...
                        l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,3);
                else
                    vj=0;
                    dvj(1)=0;
                    dvj(2)=0;
                    dvj(3)=0;
                end
                
                dphizy=(dl(2,2)*dl(1,3)-dl(1,2)*dl(2,3))*vj*wk+...
                    (l(2)*dl(1,3)-l(1)*dl(2,3))*dvj(2)*wk+...
                    (l(2)*dl(1,3)-l(1)*dl(2,3))*vj*dwk(2);
                
                dphiyz=(dl(2,3)*dl(1,2)-dl(1,3)*dl(2,2))*vj*wk+...
                    (l(2)*dl(1,2)-l(1)*dl(2,2))*dvj(3)*wk+...
                    (l(2)*dl(1,2)-l(1)*dl(2,2))*vj*dwk(3);
                
                dphixz=(dl(2,3)*dl(1,1)-dl(1,3)*dl(2,1))*vj*wk+...
                    (l(2)*dl(1,1)-l(1)*dl(2,1))*dvj(3)*wk+...
                    (l(2)*dl(1,1)-l(1)*dl(2,1))*vj*dwk(3);
                
                dphizx=(dl(2,1)*dl(1,3)-dl(1,1)*dl(2,3))*vj*wk+...
                    (l(2)*dl(1,3)-l(1)*dl(2,3))*dvj(1)*wk+...
                    (l(2)*dl(1,3)-l(1)*dl(2,3))*vj*dwk(1);
                
                dphiyx=(dl(2,1)*dl(1,2)-dl(1,1)*dl(2,2))*vj*wk+...
                    (l(2)*dl(1,2)-l(1)*dl(2,2))*dvj(1)*wk+...
                    (l(2)*dl(1,2)-l(1)*dl(2,2))*vj*dwk(1);
                
                dphixy=(dl(2,2)*dl(1,1)-dl(1,2)*dl(2,1))*vj*wk+...
                    (l(2)*dl(1,1)-l(1)*dl(2,1))*dvj(2)*wk+...
                    (l(2)*dl(1,1)-l(1)*dl(2,1))*vj*dwk(2);
 
                ph(nbas,1:3)=(dphizy-dphiyz)*asxi(1:3)+...
                    (dphixz-dphizx)*aseta(1:3)+ (dphiyx-dphixy)*aszeta(1:3);
            end
        end
    end
end
      
      if order~=0
      if nbas~=(order+1)*(order+2)*(order+3)/2
          error(message('wrong number of basis functions curl'));
      end
      end
%write(6,*)'completed curl' 