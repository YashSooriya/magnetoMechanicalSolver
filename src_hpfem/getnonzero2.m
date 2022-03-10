function [nz2,help2,irn,icn,maprow]= getnonzero2(help1,help2,nef,harwell)



display('begin reorder')
% reorder help2 
for i=1:nef
    for j=1:help1(i)-1
        for k=j+1:help1(i)
            if help2(i,j)> help2(i,k)
                store=help2(i,j);
                help2(i,j)=help2(i,k);
                help2(i,k)=store;
            end
        end
    end
end

% compute actual number of non-zero enteries
nz2=0;
for i=1:nef
    maprow(i,1)=nz2+1;
    for j=1:help1(i)
        if harwell==1
            % lower triangle for harwell
            if help2(i,j)<=i
                nz2=nz2+1;
                irn(nz2)=i;
                icn(nz2)=help2(i,j);
            end
        else
            % upper triangle for Pardiso
            if help2(i,j)>=i
                nz2=nz2+1;
                irn(nz2)=i;
                icn(nz2)=help2(i,j);
            end
        end
    end
    maprow(i,2)=nz2;
    if maprow(i,1)>maprow(i,2)
        error(message(['o enteries found for row',num2str(i)]));
    end
end

display(['The full matrix contain',num2str(nef^2),'enteries']);
display(['There are',num2str(nz2),'non-zero enteries']);