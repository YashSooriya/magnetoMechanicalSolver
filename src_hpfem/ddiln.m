function ddiln=ddiln(n,x)

ddiln=(ddln(n,x)-ddln(n-2,x))/(2*n-1);