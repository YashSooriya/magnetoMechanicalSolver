% Function to compute elemental boundary/interface integrals

function [K_UA,K_UA2,K_UA3,bhelp2,bhelp3,bhelp4,R_U,StressIntegral]= BoundTerms2(K_UA,K_UA2,K_UA3,R_U,esize,esizeH1,intfw,mue,A5,matc,unknown,order,mesh,X,faceNumber,...
    contador,area,Nhf,nm,pp,bhelp2,bhelp3,bhelp4,constantB,A2,Basis3Df,ConstantD,sigmae,probstatic,StressIntegral,jCond,ecph,eltypeCond)

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
            curlADC=ecph'*A;
                    
                
                AuxB=ecph*curlADC;
                
                % Define constant for vectorized definition
                constantB(1:3:end,1)=AuxB;
                constantB(2:3:end,2)=AuxB;
                constantB(3:3:end,3)=AuxB;
                
                % Define the components that we need
                A2(1:3:end,1)=ecph(:,1);
                A2(1:3:end,2)=ecph(:,2);
                A2(1:3:end,3)=ecph(:,3);
                A2(2:3:end,4)=ecph(:,1);
                A2(2:3:end,5)=ecph(:,2);
                A2(2:3:end,6)=ecph(:,3);
                A2(3:3:end,7)=ecph(:,1);
                A2(3:3:end,8)=ecph(:,2);
                A2(3:3:end,9)=ecph(:,3);
                
                
                % Define matrices for vectorised calculus of coupling
                % blocks
                B2=[curlADC(1) 0 0; 0 curlADC(1) 0; 0 0 curlADC(1); curlADC(2) 0 0; 0 curlADC(2) 0; 0 0 curlADC(2); curlADC(3) 0 0; 0 curlADC(3) 0; 0 0 curlADC(3)];
                B3=[curlADC(1) curlADC(2) curlADC(3); 0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3);0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3)];
                
                % Define linearised Maxwell stress tensor
                S=(1/mue)*(A2*(B2+B3)-constantB);
                % Multiply by normal vector
                constantD=S*nm';
                % Rearrange for vectorised calculus
                constantD2=reshape(constantD,3,esize);
                constantD3=repmat(constantD2,esizeH1,1);
                v = ceil([1:(esizeH1*3)]./3);
                helpMatrix=diag(Nhf(v),0);
                
                % Terms for the residual boundary term
                Constant=0.5*(curlADC'*curlADC);
                
                % Define full Maxwell stress tensor (necessary components)
                sigmae(1)=curlADC(1)*curlADC(1)-Constant;
                sigmae(2)=curlADC(2)*curlADC(2)-Constant;
                sigmae(3)=curlADC(3)*curlADC(3)-Constant;
                sigmae(4)=curlADC(1)*curlADC(2);
                sigmae(5)=curlADC(1)*curlADC(3);
                sigmae(6)=curlADC(2)*curlADC(3);
                sigmae=(1/(mue))*sigmae;
                
                % Multiply by normal vector
                ConstantD(1)=sigmae(1)*nm(1)+sigmae(4)*nm(2)+sigmae(5)*nm(3);
                ConstantD(2)=sigmae(4)*nm(1)+sigmae(2)*nm(2)+sigmae(6)*nm(3);
                ConstantD(3)=sigmae(5)*nm(1)+sigmae(6)*nm(2)+sigmae(3)*nm(3);
                
            end
        end
    end
end

% Add coupling contributions to stiffness matrix
    if contador==2
        
        K_UA2=K_UA2-(1/1.256637061435917e-06)*helpMatrix*constantD3*area*intfw(pp);
        
    elseif contador==1
        
        K_UA=K_UA-(1/1.256637061435917e-06)*helpMatrix*constantD3*area*intfw(pp);
        
    elseif contador==3
        
        K_UA3=K_UA3-(1/1.256637061435917e-06)*helpMatrix*constantD3*area*intfw(pp);
        
    end

% Add coupling contributions to residual
if probstatic==1
    %     R_U=R_U-Basis3Df*ConstantD*intfw(pp)*area;
    R_U=R_U-Basis3Df*ConstantD*intfw(pp)*area;
    StressIntegral=StressIntegral+ConstantD*intfw(pp)*area;
end


end



