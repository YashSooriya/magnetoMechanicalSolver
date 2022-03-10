function ph=basisl2(xi,eta,zeta,order,eltype)

% order = 0,1,2,3,.......

esizel2 = (order+1)*(order+2)*(order+3)/6;
% if order >=1
%     esizel2=esizel2+order*(order+1)*(order+2)/6+order*(order+1)/2
%     %(order-2)*(order-1)*order/6+(order-2)*(order-1)/2+(order-2);
% end
ph = zeros(esizel2,1);
ph(1,1) = 1;
nbas=1;

if order>=1
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
    if eltype==1
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
    
    % t values
    t(1)=l(2)+l(3);
    t(2)=l(3)+l(1);
    t(3)=l(2)+l(1);
    t(4)=l(4)+l(1);
    t(5)=l(4)+l(2);
    t(6)=l(4)+l(3);
    
    % compute derivatives of the parameterisations
    dt(1,1:3)=dl(3,1:3)+dl(2,1:3);
    dt(2,1:3)=dl(3,1:3)+dl(1,1:3);
    dt(3,1:3)=dl(2,1:3)+dl(1,1:3);
    dt(4,1:3)=dl(4,1:3)+dl(1,1:3);
    dt(5,1:3)=dl(4,1:3)+dl(2,1:3);
    dt(6,1:3)=dl(4,1:3)+dl(3,1:3);
    
    order  = order +2;
    
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
        
        for i=0:order-3
            for j=0:order-3
                for k=0:order-3
                    if i+j+k<=order-3
                        nbas=nbas+1;
                        if abs(ti(1))>1e-12
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
                        else
                            ui=0;
                            dui(1)=0;
                            dui(2)=0;
                            dui(3)=0;
                        end
                        
                        
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
                        
                        wk=l(4)*ln(k,si(3));
                        dwk(1)=dl(4,1)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,1);
                        dwk(2)=dl(4,2)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,2);
                        dwk(3)=dl(4,3)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,3);
                        
                        
                        ph(nbas,1)=dot(dwk,cross(dui,dvj));
                        
                    end
                end
            end
        end
        
        fai(1)=l(2)*dl(1,1)-dl(2,1)*l(1);
        fai(2)=l(2)*dl(1,2)-dl(2,2)*l(1);
        fai(3)=l(2)*dl(1,3)-dl(2,3)*l(1);
        
        fei = 2*cross(dl(2,:),dl(1,:));
        
        for j=0:order-3
            for k=0:order-3
                if j+k<=order-3
                    nbas=nbas+1;
                    
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
                    
                    wk=l(4)*ln(k,si(3));
                    dwk(1)=dl(4,1)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,1);
                    dwk(2)=dl(4,2)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,2);
                    dwk(3)=dl(4,3)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,3);
                    
                    ph(nbas,1)= dot(dwk,cross(dvj,fai)) +vj*dot(dwk,fei);
                    
                end
            end
        end
        
        drt = dot(dl(1,:),cross(dl(2,:),dl(3,:))) + dot(dl(2,:),cross(dl(3,:),dl(1,:)))+...
            dot(dl(3,:),cross(dl(1,:),dl(2,:)));
        
        drv = l(1)*cross(dl(2,:),dl(3,:))+ l(2)*cross(dl(3,:),dl(1,:)) + ...
            l(3)*cross(dl(1,:),dl(2,:));
        
        for k=0:order-3
            if k<=order-3
                nbas=nbas+1;
                
                wk=l(4)*ln(k,si(3));
                dwk(1)=dl(4,1)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,1);
                dwk(2)=dl(4,2)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,2);
                dwk(3)=dl(4,3)*ln(k,si(3))+l(4)*dln(k,si(3))*dsi(3,3);
                
                ph(nbas,1)= drt*wk+dot(drv,dwk);
                
            end
        end
        
    end
    
end

if nbas~=esizel2
    error(message('wrong number of basis functions'));
end