x = 1:100;

out = filter(x)

function out=filter(x)
    out=zeros(size(x));
    N = 5;
    Nofreqs = length(x);
    for i=1:Nofreqs
        if i <= N | i > Nofreqs-N
            out(i)=x(i);
        else
            out(i)=(1/(2*N +1))*sum(x(i-N:i+N));
        end
    end
end
