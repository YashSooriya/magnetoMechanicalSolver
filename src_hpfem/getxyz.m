function [x,y,z]=getxyz(xy,xi,eta,zeta)

l = zeros(1,4);
% area coordinates
l(1) = 0.5 - 0.5*xi - sqrt(3)/6*eta-sqrt(6)/12*zeta;
l(2) = 0.5 + 0.5*xi - sqrt(3)/6*eta-sqrt(6)/12*zeta;
l(3) = sqrt(3)/3*eta - sqrt(6)/12*zeta;
l(4) = sqrt(3)/(2*sqrt(2))*zeta;

x = l*xy(:,1);
y = l*xy(:,2);
z = l*xy(:,3);

