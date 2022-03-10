function ln=ln(n,x)

if x > 1
    x=1;
end

if x<-1
    x=-1;
end

ln = 0;
if n ~=0 && n ~=1
    for l=0:n
        a = dbfact(n)/(dbfact(l)*dbfact(n-l));
        b = dbfact(n)/(dbfact(n-l)*dbfact(l));
        
        if l==0
            c=1;
        else
            c= (x-1)^l;
        end
        
        if n-l ==0
            d=1;
        else
            d=(x+1)^(n-l);
        end
        ln = ln+((1/(2^n))*a*b*c*d);
    end
elseif n==1
    ln =x;
else
    ln=1;
end

