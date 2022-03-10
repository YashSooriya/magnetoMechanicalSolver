function [nmat,muz,epz,omega,mu,epl,sigma,jsrc,delta]=matinfo(matdata)


nmat = str2num(matdata{2});

muz = str2num(matdata{4});
epz = str2num(matdata{6});
omega = str2num(matdata{8});

for i = 1:nmat
    mu(i) = str2num(matdata{8+4*i-2});
    epl(i) = str2num(matdata{8+4*i-1});
    sigma(i) = str2num(matdata{8+4*i});
    jsrc(i,1:3) = str2num(matdata{8+4*nmat+2*i});
end

delta = str2num(matdata{8+6*nmat+2});
