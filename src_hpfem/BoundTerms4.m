% Function to compute elemental boundary/interface integrals

function [A,bhelp2,bhelp3,bhelp4]= BoundTerms4(esize,A5,matc,unknown,order,mesh,X,faceNumber,...
    contador,bhelp2,bhelp3,bhelp4,jCond,eltypeCond)

% Extract data from structures
globfa=mesh.face.globfa;
nelem=mesh.Nelements;
lunkv=unknown.system.unknowns;
A=zeros(esize,1);
bhelpNew=zeros(esize,1);
% Identify elements in the interface
for kk=1:nelem
    if matc(kk)==0
        for ll=1:4
            if globfa(kk,ll)==faceNumber
                
                %==========================================================
                % Renumber the DOF to match the numbering in the conducting
                % region
                %==========================================================
                eltypeNonCond=mesh.eltype(kk);
                lfedgeNonCond=get_lfedge(eltypeNonCond);
                lfedgeCond=get_lfedge(eltypeCond);
                bhelpNonCond=lunkv(1:esize,kk);
                
                % Loop over edges in a face (low order)
                for mf=1:3  % Number of edges in a face
                    mc=lfedgeCond(jCond,mf); % local edge number in conductor
                    m=lfedgeNonCond(ll,mf); % local edge number in free space
                    bhelpNew(mc)=bhelpNonCond(m);
                end
                
                % Loop over high order edges
                if order>=1
                    for vv=1:order
                        for mf=1:3 % number of edges in a face
                            mc=lfedgeCond(jCond,mf);
                            m=lfedgeNonCond(ll,mf);
                            bhelpNew(6*order+mc)=bhelpNonCond(6*order+m);
                        end
                    end
                end
               nbas=6*(order+1);
            if order>=2
                % Loop over faces
                % First type
                for hh=1:(order*order-order)/2
                    nbas=nbas+1;
                    bhelpNew((jCond-1)*(order*order-order)/2+nbas)=bhelpNonCond((ll-1)*(order*order-order)/2+nbas);
                end
                nbas=6*(order+1)+4*(order*order-order)/2;
                
                % Second type
                for hh=1:(order*order-order)/2
                    nbas=nbas+1;
                    bhelpNew((jCond-1)*(order*order-order)/2+nbas)=bhelpNonCond((ll-1)*(order*order-order)/2+nbas);
                end
                nbas=6*(order+1)+8*(order*order-order)/2;
                
                % Third type
                for hh=1:order-1
                    nbas=nbas+1;
                    bhelpNew((jCond-1)*(order-1)+nbas)=bhelpNonCond((ll-1)*(order-1)+nbas);
                end
            end
            
            %==========================================================================================
                
            % Now we define A
            A(bhelpNew>0)=X(bhelpNew(bhelpNew>0));
            for iii=1:length(A)
                if A(iii)==0
                A(iii)=A5(iii);
                end
            end

                if contador==1
                    bhelp2=bhelpNew;
                elseif contador==2
                    bhelp3=bhelpNew;
                elseif contador==3
                    bhelp4=bhelpNew;
                end
                    
                
            end
        end
    end
end





