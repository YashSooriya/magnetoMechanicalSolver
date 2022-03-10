function area=surdet(eltype,xi,eta,zeta,j,xy,lec,lfc,gorder)

% betae = zeros(1,6);

if j==1

    if eltype==1
        dxi(1,1)=0;
        dxi(1,2)=0;
        dxi(2,1)=-(1/3);
        dxi(2,2)=-(1/3);
        dxi(3,1)=(2/3);
        dxi(3,2)=-(1/3);
        dxi(4,1)=-(1/3);
        dxi(4,2)=(2/3);
    else
        dxi(1,1)=0;
        dxi(1,2)=0;
        dxi(2,1)=(2/3);
        dxi(2,2)=-(1/3);
        dxi(3,1)=-(1/3);
        dxi(3,2)=-(1/3);
        dxi(4,1)=-(1/3);
        dxi(4,2)=(2/3);
    end

elseif j==2
    
    dxi(1,1)=-(1/3);
    dxi(1,2)=-(1/3);
    dxi(2,1)=0;
    dxi(2,2)=0;
    dxi(3,1)=(2/3);
    dxi(3,2)=-(1/3);
    dxi(4,1)=-(1/3);
    dxi(4,2)=(2/3);

elseif j==3
      
    dxi(1,1)=-(1/3);
    dxi(1,2)=-(1/3);
    dxi(2,1)=(2/3);
    dxi(2,2)=-(1/3);
    dxi(3,1)=0;
    dxi(3,2)=0;
    dxi(4,1)=-(1/3);
    dxi(4,2)=(2/3);
    
else
    
    if eltype==1
        dxi(1,1)=-(1/3);
        dxi(1,2)=-(1/3);
        dxi(2,1)=(2/3);
        dxi(2,2)=-(1/3);
        dxi(3,1)=-(1/3);
        dxi(3,2)=(2/3);
        dxi(4,1)=0;
        dxi(4,2)=0;
    else
        dxi(1,1)=-(1/3);
        dxi(1,2)=-(1/3);
        dxi(2,1)=-(1/3);
        dxi(2,2)=(2/3);
        dxi(3,1)=(2/3);
        dxi(3,2)=-(1/3);
        dxi(4,1)=0;
        dxi(4,2)=0;
    end

end

% get dx/dxi vectors
dxdxi1(1:3)=dxi(1:4,1)'*xy(1:4,1:3);
dxdxi2(1:3)=dxi(1:4,2)'*xy(1:4,1:3);

% for iii=1:3
%     dxdxi1(iii)=0;
%     dxdxi2(iii)=0;
%     for jj=1:4
%         dxdxi1(iii)=dxdxi1(iii)+(xy(jj,iii)*dxi(jj,1));
%         dxdxi2(iii)=dxdxi2(iii)+(xy(jj,iii)*dxi(jj,2));
%     end
% end

      
ll(1)=0.5-(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
ll(2)=0.5+(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
ll(3)=((sqrt(3)/3)*eta)-((sqrt(6)/12)*zeta);
ll(4)=((sqrt(3)/(2*sqrt(2)))*zeta);

% paramertisations
if eltype==1
    s(1)=ll(3)-ll(2);
else
    s(1)=ll(2)-ll(3);
end
s(2)=ll(3)-ll(1);
s(3)=ll(2)-ll(1);
s(4)=ll(4)-ll(1);
s(5)=ll(4)-ll(2);
s(6)=ll(4)-ll(3);
            
% t values
t(1)=ll(2)+ll(3);
t(2)=ll(3)+ll(1);
t(3)=ll(2)+ll(1);
t(4)=ll(4)+ll(1);
t(5)=ll(4)+ll(2);
t(6)=ll(4)+ll(3);

xie = zeros(6,1);
betae = zeros(1,6);
ds = zeros(6,2);
dbetae = zeros(6,2);
      
% add edge correction
if j==1
    % face 1, we need only edges 1, 5, 6
    % face 1,edge 1
    if eltype==1
        % this is for s
        ds(1,1)=dxi(3,1)-dxi(2,1);
        ds(1,2)=dxi(3,2)-dxi(2,2);
        % this is for t
        dt(1,1)=dxi(3,1)+dxi(2,1);
        dt(1,2)=dxi(3,2)+dxi(2,2);
    else
        % this is for s
        ds(1,1)=dxi(2,1)-dxi(3,1);
        ds(1,2)=dxi(2,2)-dxi(3,2);
        % this is for t
        dt(1,1)=dxi(2,1)+dxi(3,1);
        dt(1,2)=dxi(2,2)+dxi(3,2);
    end
    betae(1)=1;
    % face 1,edge 5
    betae(5)=1;
    
    % this is for s
    ds(5,1:2)=dxi(4,1:2)-dxi(2,1:2);
    % this is for t
    dt(5,1:2)=dxi(4,1:2)+dxi(2,1:2);
    
    % face 1,edge 6
    betae(6)=1;
    
    ds(6,1:2)=dxi(4,1:2)-dxi(3,1:2);
    % this is for t
    dt(6,1:2)=dxi(4,1:2)+dxi(3,1:2);
    
    
elseif j==2
    % face 2 we need only edges 2,4,6
    % face2 edge 2
    betae(2)=1;
    
    % this is for s
    ds(2,1:2)=dxi(3,1:2)-dxi(1,1:2);
    % this is for t
    dt(2,1:2)=dxi(3,1:2)+dxi(1,1:2);
    
    % face2 edge 4
    betae(4)=1;
    
    % this is for s
    ds(4,1:2)=dxi(4,1:2)-dxi(1,1:2);
    % this is for t
    dt(4,1:2)=dxi(4,1:2)+dxi(1,1:2);
    
    % face2 edge 6
    betae(6)=1;
    
    % this is for s
    ds(6,1:2)=dxi(4,1:2)-dxi(3,1:2);
    % this is for t
    dt(6,1:2)=dxi(4,1:2)+dxi(3,1:2);
    
elseif j==3
    % face 2 we need only edges 3,4,5
    betae(3)=1;
    
    ds(3,1:2)=dxi(2,1:2)-dxi(1,1:2);
    dt(3,1:2)=dxi(2,1:2)+dxi(1,1:2);
    
    betae(4)=1;
    
    ds(4,1:2)=dxi(4,1:2)-dxi(1,1:2);
    dt(4,1:2)=dxi(4,1:2)+dxi(1,1:2);
    
    betae(5)=1;
    
    ds(5,1:2)=dxi(4,1:2)-dxi(2,1:2);
    dt(5,1:2)=dxi(4,1:2)+dxi(2,1:2);
    
else
    % face 4 we need only edges 1,2,3
    betae(1)=1;

    if eltype==1
        ds(1,1:2)=dxi(3,1:2)-dxi(2,1:2);
        dt(1,1:2)=dxi(3,1:2)+dxi(2,1:2);
    else
        ds(1,1:2)=dxi(2,1:2)-dxi(3,1:2);
        dt(1,1:2)=dxi(3,1:2)+dxi(2,1:2);
    end

    betae(2)=1;
    
    ds(2,1:2)=dxi(3,1:2)-dxi(1,1:2);
    dt(2,1:2)=dxi(3,1:2)+dxi(1,1:2);
    
    betae(3)=1;
    
    ds(3,1:2)=dxi(2,1:2)-dxi(1,1:2);
    dt(3,1:2)=dxi(2,1:2)+dxi(1,1:2);
    
end

for jj=1:6
    for k=1:gorder-1
        for iii=1:3
            if betae(jj)>0.5
                % add contributeion from this edge
                dxdxi1(iii)=dxdxi1(iii)+(-1.d0*ln(k-1,s(jj)/t(jj))*(t(jj)^k)*dt(jj,1)...
                    +ln(k,s(jj)/t(jj))*(t(jj)^k)*ds(jj,1))*lec(jj,k,iii);
                
                dxdxi2(iii)=dxdxi2(iii)+(-1.d0*ln(k-1,s(jj)/t(jj))*(t(jj)^k)*dt(jj,2)...
                    +ln(k,s(jj)/t(jj))*(t(jj)^k)*ds(jj,2) )*lec(jj,k,iii);
            end
        end
    end
end

% now add face corrections.....
betaf = zeros(4,1);
dbetaf = zeros(4,2);
xi1f = zeros(4,1);
xi2f = zeros(4,1);

if eltype==1
    sfu(1)=ll(2)-ll(3);
    sfu(4)=ll(1)-ll(2);
else
    sfu(1)=ll(3)-ll(2);
    sfu(4)=ll(1)-ll(3);
end
sfu(2)=ll(1)-ll(3);
sfu(3)=ll(1)-ll(2);

if eltype==1
    dsfu(1,1:2)=dxi(2,1:2)-dxi(3,1:2);
    dsfu(4,1:2)=dxi(1,1:2)-dxi(2,1:2);
else
    dsfu(1,1:2)=dxi(3,1:2)-dxi(2,1:2);
    dsfu(4,1:2)=dxi(1,1:2)-dxi(3,1:2);
end
dsfu(2,1:2)=dxi(1,1:2)-dxi(3,1:2);
dsfu(3,1:2)=dxi(1,1:2)-dxi(2,1:2);

      
% define tfu(m), m=1,2,3,4
if eltype==1
    tfu(1)=ll(2)+ll(3);
    tfu(4)=ll(1)+ll(2);
else
    tfu(1)=ll(3)+ll(2);
    tfu(4)=ll(1)+ll(3);
end
tfu(2)=ll(1)+ll(3);
tfu(3)=ll(1)+ll(2);

% define lf(m),m=1,2,3,4
lf(1)=ll(4);
lf(2)=ll(4);
lf(3)=ll(4);
if eltype==1
    lf(4)=ll(3);
else
    lf(4)=ll(2);
end
      
% define dlf(m,j),m=1,2,3,4,j=1,2,3
dlf(1,1:2)=dxi(4,1:2);
dlf(2,1:2)=dxi(4,1:2);
dlf(3,1:2)=dxi(4,1:2);
if eltype==1
    dlf(4,1:2)=dxi(3,1:2);
else
    dlf(4,1:2)=dxi(2,1:2);
end
    
% define sfv(m), m=1,2,3,4
if eltype==1
    sfv(1)=ll(4)-ll(2)-ll(3);
    sfv(4)=ll(3)-ll(2)-ll(1);
else
    sfv(1)=ll(4)-ll(2)-ll(3);
    sfv(4)=ll(2)-ll(3)-ll(1);
end
sfv(2)=ll(4)-ll(3)-ll(1);
sfv(3)=ll(4)-ll(2)-ll(1);
      
% define dsfv(m,j), m=1,2,3,4,j=1,2,3
if eltype==1
    dsfv(1,1:2)=dxi(4,1:2)-dxi(2,1:2)-dxi(3,1:2);
    dsfv(4,1:2)=dxi(3,1:2)-dxi(2,1:2)-dxi(1,1:2);
else
    dsfv(1,1:2)=dxi(4,1:2)-dxi(2,1:2)-dxi(3,1:2);
    dsfv(4,1:2)=dxi(2,1:2)-dxi(3,1:2)-dxi(1,1:2);
end
dsfv(2,1:2)=dxi(4,1:2)-dxi(3,1:2)-dxi(1,1:2);
dsfv(3,1:2)=dxi(4,1:2)-dxi(2,1:2)-dxi(1,1:2);
   
      
% define tfu(m), m=1,2,3,4
if eltype==1
    tfv(1)=ll(2)+ll(3)+ll(4);
    tfv(4)=ll(1)+ll(2)+ll(3);
else
    tfv(1)=ll(3)+ll(2)+ll(4);
    tfv(4)=ll(1)+ll(3)+ll(2);
end
tfv(2)=ll(1)+ll(3)+ll(4);
tfv(3)=ll(1)+ll(2)+ll(4);
      
% define dtfu(m,j), m=1,2,3,4,j=1,2,3

if eltype==1
    dtfv(1,1:2)=dxi(2,1:2)+dxi(3,1:2)+dxi(4,1:2);
    dtfv(4,1:2)=dxi(1,1:2)+dxi(2,1:2)+dxi(3,1:2);
else
    dtfv(1,1:2)=dxi(3,1:2)+dxi(2,1:2)+dxi(4,1:2);
    dtfv(4,1:2)=dxi(1,1:2)+dxi(3,1:2)+dxi(2,1:2);
end
dtfv(2,1:2)=dxi(1,1:2)+dxi(3,1:2)+dxi(4,1:2);
dtfv(3,1:2)=dxi(1,1:2)+dxi(2,1:2)+dxi(4,1:2);


if eltype==1
    dtfu(1,1:2)=dxi(2,1:2)+dxi(3,1:2);
    dtfu(4,1:2)=dxi(1,1:2)+dxi(2,1:2);
else
    dtfu(1,1:2)=dxi(3,1:2)+dxi(2,1:2);
    dtfu(4,1:2)=dxi(1,1:2)+dxi(3,1:2);
end
dtfu(2,1:2)=dxi(1,1:2)+dxi(3,1:2);
dtfu(3,1:2)=dxi(1,1:2)+dxi(2,1:2);

% we only need to consider face j!!!
kkk=0;
for k=0:gorder-3
    for kk=0:gorder-3
        if k+kk<=gorder-3
            kkk=kkk+1;
            for iii=1:3
                
                % define ui
                if abs(tfu(j))>1e-12
                    ui=iln(k+2,sfu(j)/tfu(j))*(tfu(j)^(k+2));
                else
                    ui=0;
                end
                % define vi
                if abs(tfv(j))>1e-12
                    vj=lf(j)*ln(kk,sfv(j)/tfv(j))*(tfv(j)^(kk));
                else
                    vj=0;
                end
                
                %define dui but with rpt to collapsed ccordinates
                if abs(tfu(j))>1e-12
                    dui(1)=-ln(k,sfu(j)/tfu(j))*(tfu(j)^(k+1))*dtfu(j,1)+...
                            ln(k+1,sfu(j)/tfu(j))*(tfu(j)^(k+1))*dsfu(j,1);
                    dui(2)=-ln(k,sfu(j)/tfu(j))*(tfu(j)^(k+1))*dtfu(j,2)+...
                            ln(k+1,sfu(j)/tfu(j))*(tfu(j)^(k+1))*dsfu(j,2);
                else
                    dui(1)=0;
                    dui(2)=0;
                end
                
                % define dvj but
                if abs(tfv(j))>1e-12
                    dvj(1)=dlf(j,1)*ln(kk,sfv(j)/tfv(j))*(tfv(j)^(kk))+...
                           lf(j)*dln(kk,sfv(j)/tfv(j))*((tfv(j)*dsfv(j,1)-...
                           sfv(j)*dtfv(j,1))/(tfv(j)^2))*(tfv(j)^kk)+lf(j)*...
                           ln(kk,sfv(j)/tfv(j))*kk*(tfv(j)^(kk-1))*dtfv(j,1);
                    
                    dvj(2)=dlf(j,2)*ln(kk,sfv(j)/tfv(j))*(tfv(j)^(kk))+...
                           lf(j)*dln(kk,sfv(j)/tfv(j))*((tfv(j)*dsfv(j,2)-...
                           sfv(j)*dtfv(j,2))/(tfv(j)^2))*(tfv(j)^kk)+lf(j)*...
                           ln(kk,sfv(j)/tfv(j))*kk*(tfv(j)^(kk-1))*dtfv(j,2);
                else
                    dvj(1)=0;
                    dvj(2)=0;
                end
                
                dxdxi1(iii)=dxdxi1(iii)+((dui(1)*vj)+(dvj(1)*ui))*lfc(j,kkk,iii);
                dxdxi2(iii)=dxdxi2(iii)+((dui(2)*vj)+(dvj(2)*ui))*lfc(j,kkk,iii);
                
            end
        end
    end
end


area=sqrt( ((dxdxi1(2)*dxdxi2(3))-(dxdxi1(3)*dxdxi2(2)))^2 +...
            (-1.*((dxdxi1(1)*dxdxi2(3))-(dxdxi1(3)*dxdxi2(1))))^2+...
            ((dxdxi1(1)*dxdxi2(2))-(dxdxi1(2)*dxdxi2(1)))^2 );

area=0.5*area*3;
