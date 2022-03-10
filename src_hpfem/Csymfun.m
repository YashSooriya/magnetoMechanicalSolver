function Csym = Csymfun(i,j,k,l,lambda,G)

Csym=lambda*eq(i,k)*eq(j,l)+G*(eq(i,j)*eq(k,l)+eq(i,l)*eq(k,j));

end