function ddln=ddln(n,x)

if n<=1
    ddln=0;
else
    u=((-2*n*n*x*ln(n,x))+(2*n*n*ln(n-1,x)));
    v=(2*n*(1-(x^2)));
    du=-2*n*n*ln(n,x)-2*n*n*x*dln(n,x)+2*n*n*dln(n-1,x);
    dv=-2*n*2*x;
    if v*du-u*dv==0 && v^2==0
        ddln=0;
    elseif v*du-u*dv~=0 && v^2==0
        error(message('error in gettting dln'));
    else
        ddln=(v*du-u*dv)/(v^2);
    end
end