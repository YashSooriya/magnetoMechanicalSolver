% Function to compute elemental boundary/interface integrals

function [K_UA,K_UA2,K_UA3,bhelp2,bhelp3,bhelp4]= stiffLinearBoun(K_UA,K_UA2,K_UA3,esize,esizeH1,intfw,nipf,gorder,ecphdfx,ecphdfy,ecphdfz,mue,A,matc,subFlag,unknown,order,orderH1,mesh,X,DirichletFaceBasis,faceNumber,...
    gesizet,contador,gph,area,Nhf,nm,pp,bhelp2,bhelp3,bhelp4,constantB,A2)





globfa=mesh.face.globfa;
nelem=mesh.Nelements;
intma=mesh.intma;
coord=mesh.Coordinates;




S=zeros(3*esize,3);


  for kk=1:nelem
                if matc(kk)==0
                    for ll=1:4
                        if globfa(kk,ll)==faceNumber
                             xy2 = coord(intma(kk,1:4),1:3);
                             flag2=0 ;

    mycoord2 = zeros(gesizet,3);
%    transfer coefficents to locations vertices
    mycoord2(1:4,1:3) = xy2(1:4,1:3);
    lec2=zeros(6,gorder,3);
    lfc2=zeros(4,gorder*(gorder-1)/2,3);
    if gorder > 0
        for jjj=1:6
            for p=1:gorder
                for k=1:3
                    lec2(jjj,p,k)=edgecof(glob(iii,jjj),((p-1)*3)+k);
                    if abs(edgecof(glob(iii,jjj),((p-1)*3)+k))>0.00000001
                        flag2=1;
                    end
                    mycoord2(4+jjj+6*(p-1),k)=lec2(jjj,p,k);
                end
            end
        end
        
        for jjj=1:4
            for p=1:gorder*(gorder-1)/2
                for k=1:3
                    lfc2(jjj,p,k)=facecof(globfa(iii,jjj),((p-1)*3)+k);
                    mycoord2(4+6*gorder+(jjj-1)*gorder*(gorder-1)/2+p,k)= lfc2(jjj,p,k);
                end
            end
        end
    end
    if contador==2
        
        
          bhelp3=rowfun(unknown,order,orderH1,mesh,kk,esize,esizeH1,subFlag(kk));
       
             for ii=1:esize
        row=bhelp3(ii);
        if row>0

               A(ii,1)=X(row,1);

        elseif row<0
      
               A(ii,1)=X(abs(row)+nunktEM,1);

        end
             end
        elseif contador==1
            
       bhelp2=rowfun(unknown,order,orderH1,mesh,kk,esize,esizeH1,subFlag(kk));
       
             for ii=1:esize
        row=bhelp2(ii);
        if row>0

               A(ii,1)=X(row,1);

        elseif row<0
      
               A(ii,1)=X(abs(row)+nunktEM,1);

        end
             
            
             end
elseif contador==3

            bhelp4=rowfun(unknown,order,orderH1,mesh,kk,esize,esizeH1,subFlag(kk));
             for ii=1:esize
        row=bhelp4(ii);
        if row>0

               A(ii,1)=X(row,1);

        elseif row<0
      
               A(ii,1)=X(abs(row)+nunktEM,1);

        end
             
            
             end

             end
             if mesh.eltype(kk)==1
                 ecphdfx=DirichletFaceBasis.ecphdfx1;
ecphdfy=DirichletFaceBasis.ecphdfy1;
ecphdfz=DirichletFaceBasis.ecphdfz1;
             elseif mesh.eltype(kk)==2
                 ecphdfx=DirichletFaceBasis.ecphdfx2;
ecphdfy=DirichletFaceBasis.ecphdfy2;
ecphdfz=DirichletFaceBasis.ecphdfz2;
             end
      
            
            [~,~,~,asxi2,aseta2,aszeta2,~]=jacobian_pre(flag2,gesizet,gph,mycoord2) ;
            
                                 ph(1:esize,1:3)=(ecphdfx((ll-1)*nipf+pp,1:esize)'*asxi2(1:3))+...
                (ecphdfy((ll-1)*nipf+pp,1:esize)'*aseta2(1:3))+...
                (ecphdfz((ll-1)*nipf+pp,1:esize)'*aszeta2(1:3));
            
           curlADC=ph'*A;
           AuxB=ph*curlADC;

           constantB(1:3:end,1)=AuxB;
           constantB(2:3:end,2)=AuxB;
           constantB(3:3:end,3)=AuxB;
           
           
           A2(1:3:end,1)=ph(:,1);
           A2(1:3:end,2)=ph(:,2);
           A2(1:3:end,3)=ph(:,3);
           A2(2:3:end,4)=ph(:,1);
           A2(2:3:end,5)=ph(:,2);
           A2(2:3:end,6)=ph(:,3);
           A2(3:3:end,7)=ph(:,1);
           A2(3:3:end,8)=ph(:,2);
           A2(3:3:end,9)=ph(:,3);

             
           
           B2=[curlADC(1) 0 0; 0 curlADC(1) 0; 0 0 curlADC(1); curlADC(2) 0 0; 0 curlADC(2) 0; 0 0 curlADC(2); curlADC(3) 0 0; 0 curlADC(3) 0; 0 0 curlADC(3)];
           B3=[curlADC(1) curlADC(2) curlADC(3); 0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3);0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3)];
           
           S=S+(1/mue)*(A2*(B2+B3)-constantB);
           constantD=S*nm';
           constantD2=reshape(constantD,3,esize);
           constantD3=repmat(constantD2,esizeH1,1);
           v = ceil( [1:(esizeH1*3)]./3);
           V2=Nhf(v);
           helpMatrix=diag(V2,0);                  
 
                        end
                    end
                end
            end
 
if contador==2

 K_UA2=K_UA2-helpMatrix*constantD3*area*intfw(pp);
                
elseif contador==1

 K_UA=K_UA-helpMatrix*constantD3*area*intfw(pp);
                
elseif contador==3
    
K_UA3=K_UA3-helpMatrix*constantD3*area*intfw(pp);

end

          


end



