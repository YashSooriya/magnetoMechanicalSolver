function area=surdetq(eltype,xi,eta,zeta,j,mycoord)

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

% Area coordinates on the reference element
l(1)=0.5-(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(2)=0.5+(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(3)=((sqrt(3)/3)*eta)-((sqrt(6)/12)*zeta);
l(4)=((sqrt(3)/(2*sqrt(2)))*zeta);

% replace dL_i's with dxi's
for i=1:2
    gph(1,i)=dxi(1,i)*(2*l(1)-1)+l(1)*(2*dxi(1,i));
    gph(2,i)=dxi(2,i)*(2*l(2)-1)+l(2)*(2*dxi(2,i));
    gph(3,i)=dxi(3,i)*(2*l(3)-1)+l(3)*(2*dxi(3,i));
    gph(4,i)=dxi(4,i)*(2*l(4)-1)+l(4)*(2*dxi(4,i));
    
    gph(5,i)=4*dxi(1,i)*l(2)+4*l(1)*dxi(2,i);
    gph(6,i)=4*dxi(2,i)*l(3)+4*l(2)*dxi(3,i);
    gph(7,i)=4*dxi(3,i)*l(1)+4*l(3)*dxi(1,i);
    gph(8,i)=4*dxi(1,i)*l(4)+4*l(1)*dxi(4,i);
    gph(9,i)=4*dxi(2,i)*l(4)+4*l(2)*dxi(4,i);
    gph(10,i)=4*dxi(3,i)*l(4)+4*l(3)*dxi(4,i);
end


% get dx/dxi vectors
%dxdxi1(1:3)=dxi(1:4,1)'*xy(1:4,1:3);
%dxdxi2(1:3)=dxi(1:4,2)'*xy(1:4,1:3);

dxdxi1(1:3)=gph(1:10,1)'*mycoord(1:10,1:3);
dxdxi2(1:3)=gph(1:10,2)'*mycoord(1:10,1:3);



area=sqrt( ((dxdxi1(2)*dxdxi2(3))-(dxdxi1(3)*dxdxi2(2)))^2 +...
            (-1.*((dxdxi1(1)*dxdxi2(3))-(dxdxi1(3)*dxdxi2(1))))^2+...
            ((dxdxi1(1)*dxdxi2(2))-(dxdxi1(2)*dxdxi2(1)))^2 );

area=0.5*area*3;
