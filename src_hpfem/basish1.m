function ph=basish1(esize,xi,eta,zeta,order,elotype)


ph = zeros(esize,1);
% area coordinates and there derivitives
l(1)=0.5-(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(2)=0.5+(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(3)=((sqrt(3)/3)*eta)-((sqrt(6)/12)*zeta);
l(4)=((sqrt(3)/(2*sqrt(2)))*zeta);

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
if elotype ==1
    s(1)=l(3)-l(2);
else
    s(1)=l(2)-l(3);
end
s(2)=l(3)-l(1);
s(3)=l(2)-l(1);
s(4)=l(4)-l(1);
s(5)=l(4)-l(2);
s(6)=l(4)-l(3);

% compute derivatives of the parameterisations
if elotype==1
    ds(1,1:3)=dl(3,1:3)-dl(2,1:3);
else
    ds(1,1:3)=dl(2,1:3)-dl(3,1:3);
end
ds(2,1:3)=dl(3,1:3)-dl(1,1:3);
ds(3,1:3)=dl(2,1:3)-dl(1,1:3);
ds(4,1:3)=dl(4,1:3)-dl(1,1:3);
ds(5,1:3)=dl(4,1:3)-dl(2,1:3);
ds(6,1:3)=dl(4,1:3)-dl(3,1:3);

% t values
t(1)=l(2)+l(3);
t(2)=l(3)+l(1);
t(3)=l(2)+l(1);
t(4)=l(4)+l(1);
t(5)=l(4)+l(2);
t(6)=l(4)+l(3);

% compte the derivatives of the parameterisations
dt(1,1:3)=dl(3,1:3)+dl(2,1:3);
dt(2,1:3)=dl(3,1:3)+dl(1,1:3);
dt(3,1:3)=dl(2,1:3)+dl(1,1:3);
dt(4,1:3)=dl(4,1:3)+dl(1,1:3);
dt(5,1:3)=dl(4,1:3)+dl(2,1:3);
dt(6,1:3)=dl(4,1:3)+dl(3,1:3);

% p = 1 hat functions basis functions
ph(1:4,1) = l(1:4)';
nbas = 4;
% p > 1 higher order edge basis functions
if order >=2 
    for m=1:6
        for p=1:order-1
            % note tolerence added
            if abs(t(m))>1e-12
                ph(4+m+6*(p-1),1)=iln(p+1,s(m)/t(m))*(t(m)^(p+1));
            else
                ph(4+m+6*(p-1),1)=0;
            end
        end
    end
    nbas=6*(order-1)+4;
end

% face functions
if order >= 3
    % define sfu(m), m=1,2,3,4
    if elotype==1
        sfu(1)=l(2)-l(3);
        sfu(4)=l(1)-l(2);
    else
        sfu(1)=l(3)-l(2);
        sfu(4)=l(1)-l(3);
    end
    sfu(2)=l(1)-l(3);
    sfu(3)=l(1)-l(2);
    
    % define dsfu(m,j), m=1,2,3,4, j=1,2,3
    if elotype==1
        dsfu(1,1:3)=dl(2,1:3)-dl(3,1:3);
        dsfu(4,1:3)=dl(1,1:3)-dl(2,1:3);
    else
        dsfu(1,1:3)=dl(3,1:3)-dl(2,1:3);
        dsfu(4,1:3)=dl(1,1:3)-dl(3,1:3);
    end
    dsfu(2,1:3)=dl(1,1:3)-dl(3,1:3);
    dsfu(3,1:3)=dl(1,1:3)-dl(2,1:3);
    
    % define tfu(m), m=1,2,3,4
    if elotype==1
        tfu(1)=l(2)+l(3);
        tfu(4)=l(1)+l(2);
    else
        tfu(1)=l(3)+l(2);
        tfu(4)=l(1)+l(3);
    end
    tfu(2)=l(1)+l(3);
    tfu(3)=l(1)+l(2);
    
    % define dtfu(m,j), m=1,2,3,4, j=1,2,3
    if elotype == 1
        dtfu(1,1:3)=dl(2,1:3)+dl(3,1:3);
        dtfu(4,1:3)=dl(1,1:3)+dl(2,1:3);
    else
        dtfu(1,1:3)=dl(3,1:3)+dl(2,1:3);
        dtfu(4,1:3)=dl(1,1:3)+dl(3,1:3);
    end
    dtfu(2,1:3)=dl(1,1:3)+dl(3,1:3);
    dtfu(3,1:3)=dl(1,1:3)+dl(2,1:3);
    
    % define sfv(m), m=1,2,3,4
    if elotype == 1
        sfv(1)=l(4)-l(2)-l(3);
        sfv(4)=l(3)-l(2)-l(1);
    else
        sfv(1)=l(4)-l(2)-l(3);
        sfv(4)=l(2)-l(3)-l(1);
    end
    sfv(2)=l(4)-l(3)-l(1);
    sfv(3)=l(4)-l(2)-l(1);
    
    % define dsfv(m,j) m = 1,2,3,4, j=1,2,3
    if elotype==1
        dsfv(1,1:3)=dl(4,1:3)-dl(2,1:3)-dl(3,1:3);
        dsfv(4,1:3)=dl(3,1:3)-dl(2,1:3)-dl(1,1:3);
    else
        dsfv(1,1:3)=dl(4,1:3)-dl(2,1:3)-dl(3,1:3);
        dsfv(4,1:3)=dl(2,1:3)-dl(3,1:3)-dl(1,1:3);
    end
    dsfv(2,1:3)=dl(4,1:3)-dl(3,1:3)-dl(1,1:3);
    dsfv(3,1:3)=dl(4,1:3)-dl(2,1:3)-dl(1,1:3);
    
    % define tfv(m), m=1,2,3,4
    if elotype==1
        tfv(1)=l(2)+l(3)+l(4);
        tfv(4)=l(1)+l(2)+l(3);
    else
        tfv(1)=l(3)+l(2)+l(4);
        tfv(4)=l(1)+l(3)+l(2);
    end
    tfv(2)=l(1)+l(3)+l(4);
    tfv(3)=l(1)+l(2)+l(4);
    
    % define dtfv(m,j), m=1,2,3,4, j=1,2,3
    if elotype==1
        dtfv(1,1:3)=dl(2,1:3)+dl(3,1:3)+dl(4,1:3);
        dtfv(4,1:3)=dl(1,1:3)+dl(2,1:3)+dl(3,1:3);
    else
        dtfv(1,1:3)=dl(3,1:3)+dl(2,1:3)+dl(4,1:3);
        dtfv(4,1:3)=dl(1,1:3)+dl(3,1:3)+dl(2,1:3);
    end
    dtfv(2,1:3)=dl(1,1:3)+dl(3,1:3)+dl(4,1:3);
    dtfv(3,1:3)=dl(1,1:3)+dl(2,1:3)+dl(4,1:3);
    
    % define lf(m), m=1,2,3,4
    lf(1)=l(4);
    lf(2)=l(4);
    lf(3)=l(4);
    if elotype==1
        lf(4)=l(3);
    else
        lf(4)=l(2);
    end
    
    % define dfl(m,j), m=1,2,3,4, j=1,2,3
    dlf(1,1:3)=dl(4,1:3);
    dlf(2,1:3)=dl(4,1:3);
    dlf(3,1:3)=dl(4,1:3);
    if elotype==1
        dlf(4,1:3)=dl(3,1:3);
    else
        dlf(4,1:3)=dl(2,1:3);
    end
    
    % now construct the functions for each face 
    % type 1 functions (gradients)
     nbas=6*(order-1)+4;
     for m=1:4
         for i=0:order-3
             for j=0:order-3
                 if i+j<=order-3
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
                             lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,1)-sfv(m)*dtfv(m,1))...
                             /(tfv(m)^2))*(tfv(m)^j)+...
                             lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,1);
                         
                         dvj(2)=dlf(m,2)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                             lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,2)-sfv(m)*dtfv(m,2))...
                             /(tfv(m)^2))*(tfv(m)^j)+...
                             lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,2);
                         
                         dvj(3)=dlf(m,3)*ln(j,sfv(m)/tfv(m))*(tfv(m)^(j))+...
                             lf(m)*dln(j,sfv(m)/tfv(m))*((tfv(m)*dsfv(m,3)-sfv(m)*dtfv(m,3))...
                             /(tfv(m)^2))*(tfv(m)^j)+...
                             lf(m)*ln(j,sfv(m)/tfv(m))*j*(tfv(m)^(j-1))*dtfv(m,3);
                     else
                         dvj(1)=0;
                         dvj(2)=0;
                         dvj(3)=0;
                     end
                     ph(nbas,1)=ui*vj;
                 end
             end
         end
     end  
end

if order >=4
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
    
    for i=0:order-4
        for j=0:order-4
            for k=0:order-4
                if i+j+k<=order-4
                    nbas=nbas+1;
                    if abs(ti(1))>1e-12
                        ui=iln(i+2,si(1)/ti(1))*(ti(1)^(i+2));
                        dd=((dsi(1,1)*ti(1))-(si(1)*dti(1,1)))/(ti(1)^2);
                        dui(1)=diln(i+2,si(1)/ti(1))*dd*(ti(1)^(i+2))...
                            +(iln(i+2,si(1)/ti(1))*(i+2)*(ti(1)^(i+1))*dti(1,1));     
                        dd=((dsi(1,2)*ti(1))-(si(1)*dti(1,2)))/(ti(1)^2);
                        dui(2)=diln(i+2,si(1)/ti(1))*dd*(ti(1)^(i+2))...
                            +(iln(i+2,si(1)/ti(1))*(i+2)*(ti(1)^(i+1))*dti(1,2));
                        dd=((dsi(1,3)*ti(1))-(si(1)*dti(1,3)))/(ti(1)^2);
                        dui(3)=diln(i+2,si(1)/ti(1))*dd*(ti(1)^(i+2))...
                            +(iln(i+2,si(1)/ti(1))*(i+2)*(ti(1)^(i+1))*dti(1,3));
                    else
                        ui=0;
                        dui(1)=0;
                        dui(2)=0;
                        dui(3)=0;
                    end
                    
                    if abs(ti(2))>1e-12
                        vj=l(3)*ln(j,si(2)/ti(2))*(ti(2)^j);
                        dd=((dsi(2,1)*ti(2))-(si(2)*dti(2,1)))/(ti(2)^2);
                        dvj(1)=dl(3,1)*ln(j,si(2)/ti(2))*(ti(2)^j)...
                            +l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)...
                            +l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,1);
                        dd=((dsi(2,2)*ti(2))-(si(2)*dti(2,2)))/(ti(2)^2);
                        dvj(2)=dl(3,2)*ln(j,si(2)/ti(2))*(ti(2)^j)...
                            +l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)...
                            +l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,2);
                        dd=((dsi(2,3)*ti(2))-(si(2)*dti(2,3)))/(ti(2)^2);
                        dvj(3)=dl(3,3)*ln(j,si(2)/ti(2))*(ti(2)^j)...
                            +l(3)*dln(j,si(2)/ti(2))*dd*(ti(2)^j)...
                            +l(3)*ln(j,si(2)/ti(2))*(j)*(ti(2)^(j-1))*dti(2,3);
                    else
                        vj=0;
                        dvj(1)=0;
                        dvj(2)=0;
                        dvj(3)=0;
                    end
     
                    wk=l(4)*ln(k,si(3));
                    dwk(1)=dl(4,1)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,1);
                    dwk(2)=dl(4,2)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,2);
                    dwk(3)=dl(4,3)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,3);
                    
                    ph(nbas,1)=ui*vj*wk;
                end
            end
        end
    end       
end


 if nbas~=(order+3)*(order+2)*(order+1)/6
     error(message('wrong number of basis functions'));
 end