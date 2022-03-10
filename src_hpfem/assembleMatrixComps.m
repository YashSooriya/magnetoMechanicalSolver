function [I,J,Kv,nmst]=assembleMatrixComps(K,I,J,Kv,nmst,bhelpRows,bhelpCols,sizeRows,sizeCols)

%-------------------------------------------------------------------------
% Local to global matrix assembly
%-------------------------------------------------------------------------
% Assemble linear system in vector format
for j = 1:sizeRows %econt
    row = bhelpRows(j);
    if row>0
        for k = 1:sizeCols  %econt
            col = bhelpCols(k);
            if col>0
                nmst=nmst+1;
                I(nmst) = row;
                J(nmst) = col;
                Kv(nmst) =Kv(nmst)+ K(j,k);
                
            end
        end
        % Add right hand side source terms

    end
end