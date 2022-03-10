 function [b,t,w]=gaussq1(kind,n,alpha,beta,kpts,endpts)


[b, t, muzero]=class(kind, n, alpha, beta);


if kpts ~=0
    if kpts==2
        gam = solve(endpts(1), n, t, b);
        t1 = ((endpts(1) - endpts(2))/(solve(endpts(2), n, t, b) - gam));
        b(n-1) = sqrt(t1);
        t(n) = endpts(1) + gam*t1;
    else
        t(n) = solve(endpts(1), n, t, b)*b(n-1)^2 + endpts(1);
    end
end

w(1) = 1.0d0;
for i = 2:n
    w(i) = 0.0d0;
end

[b,t,w,ierr]=gausq2(n, t, b, w);

for i = 1:n
    w(i) = muzero * w(i) * w(i);
end

return
end

function sval = solve(shift, n, a, b)

alpha = a(1) - shift;
nm1 = n - 1;
for i = 2:nm1
    alpha = a(i) - shift - b(i-1)^2/alpha;
end
sval = 1.0d0/alpha;
return
end

      
function [b, a, muzero]=class(kind, n, alpha, beta)

pi = 4.0d0 * atan(1.0d0);
nm1 = n - 1;

if kind ==1
    muzero = 2.0d0;
    for i = 1:nm1
        a(i) = 0.0d0;
        abi = i;
        b(i) = abi/sqrt(4*abi*abi - 1.0d0);
    end
    a(n) = 0.0d0;
    return
elseif kind ==2
    muzero = pi;
    for i = 1:nm1
        a(i) = 0.0d0;
        b(i) = 0.5d0;
    end
    b(1) = sqrt(0.5d0);
    a(n) = 0.0d0;
    return
elseif kind==3
    muzero = pi/2.0d0;
    for i = 1:nm1
        a(i) = 0.0d0;
        b(i) = 0.5d0;
    end
    a(n) = 0.0d0;
    return
elseif kind==4
    muzero = sqrt(pi);
    for i = 1:nm1
        a(i) = 0.0d0;
        b(i) = sqrt(i/2.0d0);
    end
    a(n) = 0.0d0;
    return
elseif kind==5
    ab = alpha + beta;
    abi = 2.0d0 + ab;
    muzero = 2.0d0^(ab +1.0d0)*gamma(alpha+1.0d0)*gamma(beta+1.0d0)/gamma(abi);
    a(1) = (beta - alpha)/abi;
    b(1) = sqrt(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/((abi + 1.0d0)*abi*abi));
    a2b2 = beta*beta - alpha*alpha;
    for i = 2:nm1
        abi = 2.0d0*i + ab;
        a(i) = a2b2/((abi - 2.0d0)*abi);
        b(i) = sqrt(4.0d0*i*(i+alpha)*(i+beta)*(i+ab)/((abi*abi-1)*abi*abi));
    end
    abi = 2.0d0*n + ab;
    a(n) = a2b2/((abi - 2.0d0)*abi);
    return
elseif kind==6
    muzero = gamma(alpha + 1.0d0);
    for i = 1:nm1
        a(i) = 2.0d0*i - 1.0d0 + alpha;
        b(i) = sqrt(i*(i + alpha));
    end
    a(n) = 2.0d0*n - 1 + alpha;
    return
end
end


function [e,d,z,ierr]=gausq2(n, d, e, z) % n t b w

machep=eps;

ierr = 0;

if n==1
    return
end

e(n) = 0.0d0;
for l = 1:n % 240
    j = 0;
    % :::::::::: look for small sub-diagonal element ::::::::::
    while 1 %105
        for m = l:n % 110
            if m==n
                break
            end
            if (abs(e(m)) <= machep * (abs(d(m)) + abs(d(m+1))))
                break
            end
        end %110
        p = d(l); % 120
        if (m==l)
            break
        end
        if (j == 30)
            ierr = l;
            return
        end
        j = j + 1;
        % :::::::::: form shift ::::::::::
        g = (d(l+1) - p) / (2.0d0 * e(l));
        r = sqrt(g*g+1.0d0);
        if g >=0
            sig = 1;
        else
            sig=-1;
        end
        g = d(m) - p + e(l) / (g + abs(r)*sig); %dsign(r, g)
        s = 1.0d0;
        c = 1.0d0;
        p = 0.0d0;
        mml = m - l;
        % :::::::::: for i=m-1 step -1 until l do -- ::::::::::
        for ii = 1:mml % 200
            i = m - ii;
            f = s * e(i);
            b = c * e(i);
            if (abs(f) < abs(g))
                s = f / g;
                r = sqrt(s*s+1.0d0);
                e(i+1) = g * r;
                c = 1.0d0 / r;
                s = s * c;
            else
                c = g / f;
                r = sqrt(c*c+1.0d0);
                e(i+1) = f * r;
                s = 1.0d0 / r;
                c = c * s;
            end
            g = d(i+1) - p;
            r = (d(i) - g) * s + 2.0d0 * c * b;
            p = s * r;
            d(i+1) = g + p;
            g = c * r - b;
            % :::::::::: form first component of vector ::::::::::
            f = z(i+1);
            z(i+1) = s * z(i) + c * f;
            z(i) = c * z(i) - s * f;
        end % 200
        d(l) = d(l) - p;
        e(l) = g;
        e(m) = 0.0d0;
    end %while
end % 240
% :::::::::: order eigenvalues and eigenvectors ::::::::::
for ii = 2:n % 300
    i = ii - 1;
    k = i;
    p = d(i);
    for j = ii:n % 260
        if (d(j) >= p)
            continue
        end
        k = j;
        p = d(j);
    end %260
    if (k == i)
        continue
    end
    d(k) = d(i);
    d(i) = p;
    p = z(i);
    z(i) = z(k);
    z(k) = p;
end% 300

return
% :::::::::: last card of gausq2 ::::::::::
end
