function [unknown]=elemenntUnknownNumbering(mesh,ProblemData)

nelem=mesh.Nelements;
order=ProblemData.order;
orderH1=ProblemData.orderH1;
esize=ProblemData.esizet;
esizeH1=ProblemData.esizeH1;
subFlag=mesh.subFlag;

[unknown]=nounk(mesh,ProblemData);
lunkv=zeros(esizet+3*esizeH1,nelem);

for i=1:nelem
    
    [lunkv]=rowfun(unknown,order,orderH1,mesh,i,esize,esizeH1,subFlag,lunkv);
end
    
    