% Function to compute Gauss-Legendre quadrature points on reference
% tetrahedron

function [intxi,inteta,intzeta,intw,nipv]=intpoints(nipd,maxnip)

n = nipd;

if n^3 > maxnip
    error(message('increase the dimension maxnip'));
end
%display('Integration points and weights now calculating using NETLIB routines')

% select -1 to 1 integration
kpts = 0;
endpts(1)=-1;
endpts(2)=1;

% compute weights for a face of a tetrahedra
% part 1 
% int_-1^+1 (1+x)^2 dx
alpha = 0;
beta = 2;
kind =5;
[bb,x,a]=gaussq1(kind,n,alpha,beta,kpts,endpts);
% now t = (1-x)/2 and a = a/8
t = (1-x)./2;
a = a./8;
%display('complete part 1')

% part 2
% int_-1^+1 (1+x) dx
alpha = 0;
beta = 1;
kind =5;
[bb,x,b]=gaussq1(kind,n,alpha,beta,kpts,endpts);
% now u = (1-x)/2 and b = b/4
u = (1-x)./2;
b = b./4;
%display('complete part 2')

% part 3
% int_-1^+1 dx
kind=1;
[bb,x,c]=gaussq1(kind,n,alpha,beta,kpts,endpts);
% now v = (x+1)/2 and c=c/2
v = (1+x)./2;
c = c./2;
%display('complete part 3')


% compute weights and locations
p=0;
for i=1:n
    for j=1:n
        for k=1:n
            p=p+1;
            intxi(p)=t(i);
            inteta(p)=u(j)*(1-t(i));
            intzeta(p)=v(k)*(1-t(i))*(1-u(j));
            intw(p)=a(i)*b(j)*c(k);
        end
    end
end
nipv=p;

% now alter points so that they are contained in Joe's tetrahedra
xp(1)=-1.;
xp(2)=1.;
xp(3)=0.;
xp(4)=0.;

yp(1)=0.;
yp(2)=0.;
yp(3)=sqrt(3.);
yp(4)=sqrt(3.)/3.;

zp(1)=0.;
zp(2)=0.;
zp(3)=0.;
zp(4)=2*(sqrt(2.)/sqrt(3.));

for p=1:nipv
    xi=0.;
    eta=0.;
    zeta=0.;
    l(1)=1-intxi(p)-inteta(p)-intzeta(p);
    l(2)=intxi(p);
    l(3)=inteta(p);
    l(4)=intzeta(p);
    for j=1:4
        xi=xi+(l(j)*xp(j));
        eta=eta+(l(j)*yp(j));
        zeta=zeta+(l(j)*zp(j));
    end
    intxi(p)=xi;
    inteta(p)=eta;
    intzeta(p)=zeta;
end

% alter weighting functions
for i=1:nipv
    intw(i)=intw(i)*(4*sqrt(2.));
end

% % compute weights and locations
% p = 0;
% for i=1:n
%     for j=1:n
%         for k = 1:n
%             p=p+1;
%             intxi(p)=t(i);
%             inteta(p)=u(j)*(1-t(i));
%             intzeta(p)=v(k)*(1-t(i))*(1-u(j));
%             intw(p) = a(i)*b(j)*c(k);
%         end
%     end
% end
% nipv=p;
% 
% % now alter points so that they are contained in Joe`s tetrahedra
% xp(1) = -1.;
% xp(2) = 1.;
% xp(3) = 0.;
% xp(4) = 0.;
% 
% yp(1) = 0.;
% yp(2) = 0.;
% yp(3) = sqrt(3.);
% yp(4) = sqrt(3.)/3.;
% 
% zp(1) = 0.;
% zp(2) = 0.;
% zp(3) = 0.;
% zp(4) =2*(sqrt(2.)/sqrt(3.));
% 
% for p = 1:nipv
%     xi=0.;
%     eta=0.;
%     zeta=0.;
%     l(1)=1-intxi(p)-inteta(p)-intzeta(p);
%     l(2)=intxi(p);
%     l(3)=inteta(p);
%     l(4)=intzeta(p);
%     for j=1:4
%         xi=xi+(l(j)*xp(j));
%         eta=eta+(l(j)*yp(j));
%         zeta=zeta+(l(j)*zp(j));
%     end
%     
%     intxi(p) = xi;
%     inteta(p) = eta;
%     intzeta(p) = zeta;
% end
% 
% intw(1:nipv) = intw(1:nipv)*(4*sqrt(2));
% for i = 1:nipv
%     intw(i) = intw(i)*(4*sqrt(2));
% end
%display('The computed points and weghts are ......')