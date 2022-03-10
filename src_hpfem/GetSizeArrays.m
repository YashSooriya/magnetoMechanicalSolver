function [ProblemData]=GetSizeArrays(ProblemData,mesh)

order=ProblemData.order;
orderH1=ProblemData.orderH1;


nelem=mesh.Nelements;

% assain the order to elements

    orderel = order*ones(nelem,1);


% size of local arrays exc interiors
esize = 6*(order+1);
if order >=2
    esize = esize +(12*(order-1))+(4*(order-1)*(order-2));
end


% size of local arrays including interiors
esizet = esize;
if order>=3
    esizet=esizet+2*(order-1)*(order-2)+(order-2)*(order-1)*(order-3)/2 ;
end


% %Size of local arrays including H1
esizeH1=0;
% Vertex unknowns
esizeH1=esizeH1+4;


%High order unknowns
if orderH1>=2
    esizeH1=esizeH1+4*((orderH1-1)^2-(orderH1-1))/2+(orderH1-3)*(orderH1-2)*(orderH1-1)/6+6*(orderH1-1);
end



ProblemData.esizet=esizet;
ProblemData.esize=esize;
ProblemData.esizeH1=esizeH1;
ProblemData.orderel=orderel;