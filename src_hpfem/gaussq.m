function [b,t,w]=gaussq(kind,n,alpha,beta,kpts,endpts)

% Input
% kind = an integer between 1 and 6 giving the type of quadrature rule;
% kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
% kind = 2:  chebyshev quadrature of the first kind
%            w(x) = 1/sqrt(1 - x*x) on (-1, +1)
% kind = 3:  chebyshev quadrature of the second kind
%            w(x) = sqrt(1 - x*x) on (-1, 1)
% kind = 4:  hermite quadrature, w(x) = exp(-x*x) on
%            (-infinity, +infinity)
% kind = 5:  jacobi quadrature, w(x) = (1-x)^alpha * (1+x)^beta on (-1, 1),
%            alpha, beta> -1.note: kind=2 and 3 are a special case of this.
% kind = 6:  generalized laguerre quadrature, w(x) = exp(-x)*x^alpha 
%            on (0, +infinity), alpha >-1
% n        the number of points used for the quadrature rule
% alpha    real parameter used only for gauss-jacobi and gauss-
%          laguerre quadrature (otherwise use 0.d0).
% beta     real parameter used only for gauss-jacobi quadrature--
%          (otherwise use 0.d0)
% kpts     (integer) normally 0, unless the left or right end-point (or both) 
%           of the interval is required to be anode (this is called 
%           gauss-radau or gauss-lobatto quadrature).  then kpts is 
%           the number of fixed endpoints (1 or 2).
% endpts    real array of length 2.  contains the values of
%           any fixed endpoints, if kpts = 1 or 2.
% b        real scratch array of length n

% Output
% t        will contain the desired nodes.
% w        will contain the desired weights w(j).



w = zeros(1,n);

[b,t,muzero]=class(kind,n,alpha,beta);
if kpts ~=0 && kpts~=2
    % if kpts==1, only t(n) must be changed
    t(n) = solve(endpts(1),n,t,b)*b(n-1)^2 +endpts(1);
elseif kpts == 2
    % if kpts==2, t(n) and b(n-1) must be recomputed
    gam = solve(endpts(1),n,t,b);
    t1 = ((endpts(1)-endpts(2))/(solve(endpts(2),n,t,b)-gam));
    b(n-1) = sqrt(t1);
    t(n) = endpts(1) + gam*t1;
end

w(1,1) = 1;
% w(1,2:n) = zeros(1,n-1);

[t,w,ierr]=gausq2(n,t,b,w);

w = muzero*(w.*w);

return
end


function sol = solve(shift,n,a,b)

alpha = a(1)-shift;
nm1 = n-1;
for i = 2:nm1
    alpha = a(i)-shift-(b(i-1)^2)/alpha;  
end

sol = 1/alpha;
return
end



function [b,a,muzero]=class(kind,n,alpha,beta)

nm1 = n-1; 
a = zeros(1,n);
b = zeros(1,n);
switch kind
    case 1
        muzero = 2;
        for i = 1:nm1
            a(i) = 0;
            abi = i;
            b(i) = abi/sqrt(4*abi*abi - 1);
        end
        a(n) = 0;
        return
    case 2
        muzero = pi;
        for i = 1:nm1;
            a(i) = 0;
            b(i) = 0.5;
        end
        b(1) = sqrt(0.5);
        a(n) = 0;
        return
    case 3 
        muzero = pi/2;
        for i = 1:nm1
            a(i) = 0;
            b(i) = 0.5;
        end
        a(n) = 0;
        return
    case 4
        muzero = sqrt(pi);
        for i = 1:nm1
            a(i) = 0;
            b(i) = sqrt(i/2);
        end
        a(n) = 0;
        return
    case 5
        ab = alpha+beta;
        abi = 2 +ab;
        muzero = 2^(ab+1)*gamma(alpha+1)*gamma(beta+1)/gamma(abi);
        a(1) = (beta-alpha)/abi;
        b(1) = sqrt(4*(1+alpha)*(1+beta)/((abi+1)*abi*abi));
        a2b2 = beta*beta - alpha*alpha;
        for i = 2:nm1
            abi = 2*i +ab;
            a(i) = a2b2/((abi-2)*abi);
            b(i) = sqrt(4*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-1)*abi*abi));
        end
        abi = 2*n+ab;
        a(n) = a2b2/((abi-2)*abi);
        return
    case 6
        muzero = gamma(alpha+1);
        for i = 1:nm1
            a(i) = 2*i-1+alpha;
            b(i) = sqrt(i*(i+alpha));
        end
        a(n) = 2*n - 1 +alpha;
        return
end

end
    
function [d,z,ierr]=gausq2(n,d,e,z) % n t b w

% Input
% n    is the order of the matrix;
% d    contains the diagonal elements of the input matrix;
% e    contains the subdiagonal elements of the input matrix
%      in its first n-1 positions.  e(n) is arbitrary;
% z    contains the first row of the identity matrix.
%
% Output
% d    contains the eigenvalues in ascending order.  if an
%      error exit is made, the eigenvalues are correct but
%      unordered for indices 1, 2, ..., ierr-1;
% e    has been destroyed;
% z    contains the first components of the orthonormal eigenvectors
%      of the symmetric tridiagonal matrix.  if an error exit is made,
%      z contains the eigenvectors associated with the stored eigenvalues;
% ierr is set to
%      zero for normal return,
%      j if the j-th eigenvalue has not been determined after 30 iterations.

machep = eps;

ierr = 0;

if n ==1
    return
end

e(n) = 0;
for l = 1:n %240
    j = 0;
    % look for small sub-diagonal element
    while 1
        flagm=0;
        for m = l:n % 110
            if (m ==n) || (abs(e(m))<=machep*(abs(d(m))+abs(d(m+1))))
                p = d(l);
                if m == l
                    flagm = 1;
                    break
                end
                if j == 30
                    ierr = l;
                    return
                end
                j = j+1;
                g = (d(l+1)-p)/(2*e(l));
                r = sqrt(g*g+1);
                if g>=0
                    sig = 1;
                else
                    sig = -1;
                end
                g = d(m) - p + e(l)/(g+sig*abs(r));
                s = 1;
                c = 1;
                p = 0;
                mml = m-l;
                
                for ii = 1:mml %200
                    i = m-ii;
                    f = s*e(i);
                    b = c*e(i);
                    if abs(f) < abs(g)
                        s = f/g;
                        r = sqrt(s*s+1);
                        e(i+1) = g*r;
                        c = 1/r;
                        s = s*c;
                    else
                        c = g/f;
                        r = sqrt(c*c+1);
                        e(i+1) = f*r;
                        s = 1/r;
                        c = c*s;
                    end
                    g = d(i+1) -p;
                    r = (d(i)-g)*s +2*c*b;
                    p = s*r;
                    d(i+1) = g+p;
                    g = c*r - b;
                    
                    f = z(i+1);%z(i) + c*f;
                    z(i+1) = s*z(i)+c*f;
                    z(i) = c*z(i)-s*f;
                end %200
                d(l) = d(l)-p;
                e(l) =g;
                e(m) = 0;
                break
            end
        end %110
        if flagm ==1
            break
        end
    end% while
end%240

for ii = 2:n % 300
    i = ii-1;
    k = i;
    p = d(i);
    
    for j = ii:n %260
        if d(j)>=p
            continue  
        end
        k=j;
        p=d(j);
    end%260
    if k==i
        continue
    end
    d(k)=d(i);
    d(i)=p;
    p = z(i);
    z(i) = z(k);
    z(k) = p;
end%300

return
end
    