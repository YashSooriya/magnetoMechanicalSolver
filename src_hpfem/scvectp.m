function tp = scvectp(v,w,z)

u(1) = w(2)*z(3) - w(3)*z(2);
u(2) = -1*(w(1)*z(3)-w(3)*z(1));
u(3) = w(1)*z(2) - w(2)*z(1);

tp = v*u';