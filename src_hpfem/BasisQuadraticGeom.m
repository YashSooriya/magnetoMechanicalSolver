function dN= BasisQuadraticGeom(xi,eta,zeta)

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

dN=zeros(10,3);

for i=1:3
    dN(1,i)=dl(1,i)*(2*l(1)-1)+2*l(1)*dl(1,i);
    dN(2,i)=dl(2,i)*(2*l(2)-1)+2*l(2)*dl(2,i);
    dN(3,i)=dl(3,i)*(2*l(3)-1)+2*l(3)*dl(3,i);
    dN(4,i)=dl(4,i)*(2*l(4)-1)+2*l(4)*dl(4,i);
    dN(5,i)=4*dl(1,i)*l(2)+4*l(1)*dl(2,i);
    dN(6,i)=4*dl(2,i)*l(3)+4*l(2)*dl(3,i);
    dN(7,i)=4*dl(3,i)*l(1)+4*l(3)*dl(1,i);
    dN(8,i)=4*dl(1,i)*l(4)+4*l(1)*dl(4,i);
    dN(9,i)=4*dl(2,i)*l(4)+4*l(2)*dl(4,i);
    dN(10,i)=4*dl(3,i)*l(4)+4*l(3)*dl(4,i);
end

