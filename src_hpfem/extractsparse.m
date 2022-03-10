% function to extract block and store in sparse block diagonal format

function [Xout]=extractsparse(X,unk,flag)
% input list of local numbering in the form
%unk (nentities, size of block)

[nentities,nblock]=size(unk);

I=zeros(nentities*nblock,1);
J=I;
Xz=I;
nz=0;
for i=1:nentities
    for j=1:nblock
        row=unk(i,j)-flag;
        if row>0
            for k=1:nblock
                
                col=unk(i,k)-flag;
                if col > 0
                    nz=nz+1;
                    I(nz)=row;
                    J(nz)=col;
                    Xz(nz)=X(row,col);
                end
            end
        end
    end
end
Xout=sparse(I(1:nz),J(1:nz),Xz(1:nz));
