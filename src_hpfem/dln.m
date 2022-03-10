function dln=dln(n,x)

if x>1 
    x=1;
end

if x < -1
    x=-1;
end

if n <=0
    dln=0;
else
    if 1-x^2 ==0
        dln=0;
    else
        dln=((-2.*n*n*x*ln(n,x))+(2*n*n*ln(n-1,x)))/(2*n*(1-(x^2)));
    end
end