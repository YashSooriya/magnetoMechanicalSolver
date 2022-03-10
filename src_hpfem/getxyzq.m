function [x,y,z]= getxyzq(xy,xi,eta,zeta)
	
l(1)=0.5-(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(2)=0.5+(0.5*xi)-((sqrt(3)/6)*eta)-((sqrt(6)/12)*zeta);
l(3)=((sqrt(3)/3)*eta)-((sqrt(6)/12)*zeta);
l(4)=((sqrt(3)/(2*sqrt(2)))*zeta);

n=zeros(1,10);
n(1)=l(1)*(2*l(1)-1);
n(2)=l(2)*(2*l(2)-1);
n(3)=l(3)*(2*l(3)-1);
n(4)=l(4)*(2*l(4)-1);
n(5)=4*l(1)*l(2);
n(6)=4*l(2)*l(3);
n(7)=4*l(3)*l(1);
n(8)=4*l(1)*l(4);
n(9)=4*l(2)*l(4);
n(10)=4*l(3)*l(4);

x = n*xy(1:10,1);
y = n*xy(1:10,2);
z = n*xy(1:10,3);



