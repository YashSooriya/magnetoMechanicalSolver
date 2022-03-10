% Function to compute elemental boundary/interface integrals

function [K_UA,K_UA2,K_UA3,bhelp2,bhelp3,bhelp4,R_U,StressIntegral]= BoundTerms(K_UA,K_UA2,K_UA3,R_U,esize,esizeH1,intfw,nipf,gorder,ecphdfx,ecphdfy,ecphdfz,mue,A5,matc,subFlag,unknown,order,orderH1,mesh,X,Basis,faceNumber,...
    gesizet,contador,gph,area,Nhf,nm,pp,bhelp2,bhelp3,bhelp4,constantB,A2,Basis3Df,ConstantD,sigmae,probstatic,StressIntegral,x,y,z,j,curlADC,ph2,mycoord,ecph)

% Extract data from structures
globfa=mesh.face.globfa;
nelem=mesh.Nelements;
intma=mesh.intma;
coord=mesh.Coordinates;
lunkv=unknown.system.unknowns;
pp2=0;
A=zeros(esize,1);
% Identify elements in the interface
for kk=1:nelem
    if matc(kk)==0
        for ll=1:4
            if globfa(kk,ll)==faceNumber
                xy2 = coord(intma(kk,1:4),1:3);
                flag2=0 ;
                mycoord2 = zeros(gesizet,3);
                %    transfer coefficents to locations vertices
                mycoord2(1:4,1:3) = xy2;
                if gorder > 0
                    mycoord2=mesh.mycoord(:,:,kk);
                end
                
                % Extract unknown numbering and create help arrays for new
                % numbering
                if contador==2
                    bhelp3=lunkv(1:esize,kk);
                    A(bhelp3>0)=X(bhelp3(bhelp3>0),1);
                elseif contador==1
                    bhelp2=lunkv(1:esize,kk);
                    A(bhelp2>0)=X(bhelp2(bhelp2>0),1);
                elseif contador==3
                    bhelp4=lunkv(1:esize,kk);
                    A(bhelp4>0)=X(bhelp4(bhelp4>0),1);
                end
                
                % Extract face basis
                if mesh.eltype(kk)==1
                    ecphdfx=Basis.Face.ecphdfx1;
                    ecphdfy=Basis.Face.ecphdfy1;
                    ecphdfz=Basis.Face.ecphdfz1;
                    gphfx=Basis.Face.gphfx1;
                    gphfy=Basis.Face.gphfy1;
                    gphfz=Basis.Face.gphfz1;
                    phh1f=Basis.Face.phh1f1;
                elseif mesh.eltype(kk)==2
                    ecphdfx=Basis.Face.ecphdfx2;
                    ecphdfy=Basis.Face.ecphdfy2;
                    ecphdfz=Basis.Face.ecphdfz2;
                    gphfx=Basis.Face.gphfx2;
                    gphfy=Basis.Face.gphfy2;
                    gphfz=Basis.Face.gphfz2;
                    phh1f=Basis.Face.phh1f2;
                end
                
                

   

                
                posIn=[x;y;z];
                for ppp=1:nipf 
                ph1(1:gesizet,1)=phh1f((ll-1)*nipf+ppp,1:gesizet)';
                aaa=ph1(2);
                ph1(2)=ph1(3);
                ph1(3)=aaa;
                [x2,y2,z2]= getxyzcu_pre(ph1,mycoord2,gesizet);
                posOut=[x2;y2;z2];
                if norm(posOut-posIn)<1e-8
                    pp2=ppp;
                end
                end



                 
                gph(1:gesizet,1)=gphfx((ll-1)*nipf+pp2,1:gesizet)';
                gph(1:gesizet,2)=gphfy((ll-1)*nipf+pp2,1:gesizet)';
                gph(1:gesizet,3)=gphfz((ll-1)*nipf+pp2,1:gesizet)';
                % Mapping
                [~,~,~,asxi2,aseta2,aszeta2,~]=jacobian_pre(flag2,gesizet,gph,mycoord2);
                
                ph(1:esize,1:3)=(ecphdfx((ll-1)*nipf+pp2,1:esize)'*asxi2(1:3))+(ecphdfy((ll-1)*nipf+pp2,1:esize)'*aseta2(1:3))+(ecphdfz((ll-1)*nipf+pp2,1:esize)'*aszeta2(1:3));
                
                % Compute curl(A)
                curlADC=ph'*A;
                curlADC2=[1;0;0];
                curlADC3=curlADC2-curlADC;
                if norm(curlADC3)>1e-8
                    aaaaa=1;
                end
                if j~=ll
                    aaaaa=1;
                end
                
                AuxB=ph*curlADC;
                
                % Define constant for vectorized definition
                constantB(1:3:end,1)=AuxB;
                constantB(2:3:end,2)=AuxB;
                constantB(3:3:end,3)=AuxB;
                
                % Define the components that we need
                A2(1:3:end,1)=ph(:,1);
                A2(1:3:end,2)=ph(:,2);
                A2(1:3:end,3)=ph(:,3);
                A2(2:3:end,4)=ph(:,1);
                A2(2:3:end,5)=ph(:,2);
                A2(2:3:end,6)=ph(:,3);
                A2(3:3:end,7)=ph(:,1);
                A2(3:3:end,8)=ph(:,2);
                A2(3:3:end,9)=ph(:,3);
                
                
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



