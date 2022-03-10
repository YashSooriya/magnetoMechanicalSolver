% Function for computing Gauss-Legendre points for integrating over
% triangular faces.

function [intxi,inteta,intw,nipf]=intpointf(nipd,maxnip)

%display('Computing the weights and points for the integral')
%display('int_0^1 int_0^1-x f(xy) dy dx')

n = nipd;

if n^2 > maxnip
    error(message('increase the dimension maxnip'));
end
%display('Integration points and weights now calculating using NETLIB routines')

% select -1 to 1 integration
kpts = 0;
endpts(1)=-1;
endpts(2)=1;

% compute weights for a face of a tetrahedra
% part 1
% int_-1^+1 (1+x) dx
alpha = 0;
beta = 1;
kind =5;

[bb,x,b]=gaussq1(kind,n,alpha,beta,kpts,endpts);
% now u = (1-x)/2 and b = b/4
for i=1:n
    u(i) = (1-x(i))/2;
    b(i) = b(i)/4;
end
display('complete part 1')

% part 2
% int_-1^+1 dx
kind=1;
[bb,x,c]=gaussq1(kind,n,alpha,beta,kpts,endpts);
% now v = (x+1)/2 and c=c/2
for i = 1:n
    v(i) = (1+x(i))/2;
    c(i) = c(i)/2;
end
display('complete part 2')

% compute weights and locations
p = 0;
for i=1:n
    for j=1:n
        p=p+1;
        intxi(p)=u(i);
        inteta(p)=v(j)*(1-u(i));
        intw(p) = b(i)*c(j);
    end
end
nipf=p;

for i = 1:nipf
    intw(i) = intw(i)*2;
end

