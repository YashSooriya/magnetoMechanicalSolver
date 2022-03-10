function [mesh]=ConductorBodiesFlag(probdata,mesh)

mat=mesh.mat;
nelem=mesh.Nelements;
matc=zeros(nelem,1);
nelemCond=0;

for i=1:nelem
    for j=probdata.matr.matcond
        if mat(i)==j
            nelemCond=nelemCond+1;
            matc(i)=1;
        end
    end
end
subFlag=zeros(nelem,1);
for i=1:nelem
    for j=probdata.matr.subdom_mech
        if mat(i)==j
            subFlag(i)=2;
        end
    end
end
mesh.matc=matc;
mesh.subFlag=subFlag;

disp(['number of elements in conductor = ',num2str(nelemCond)]);
disp(['Total number of elements = ',num2str(nelem)]);