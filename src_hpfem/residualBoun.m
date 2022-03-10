% Function to compute elemental boundary/interface integrals

function [R_U]= residualBoun(R_U,esize,esizeH1,intfw,nipf,gorder,ecphdfx,ecphdfy,ecphdfz,A,matc,subFlag,unknown,order,orderH1,mesh,X,DirichletFaceBasis,faceNumber,gesizet,gph,area,nm,sigmae,pp,Basis3Df,ConstantD,mue)



globfa=mesh.face.globfa;
nelem=mesh.Nelements;
intma=mesh.intma;
coord=mesh.Coordinates;

       
            
            for kk=1:nelem
                if matc(kk)==0
                    for ll=1:4
                        if globfa(kk,ll)==faceNumber
                             xy2 = coord(intma(kk,1:4),1:3);
                             flag2=0 ;
    mycoord2 = zeros(gesizet,3);
    %transfer coefficents to locations vertices
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
            
       bhelp=rowfun(unknown,order,orderH1,mesh,kk,esize,esizeH1,subFlag(kk));
       
             for ii=1:esize
        row=bhelp(ii);
        if row>0

               A(ii,1)=X(row,1);

        elseif row<0
      
               A(ii,1)=X(abs(row)+nunktEM,1);

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
            
            Constant=0.5*(curlADC'*curlADC);
            
            sigmae(1)=curlADC(1)*curlADC(1)-Constant;
            sigmae(2)=curlADC(2)*curlADC(2)-Constant;
            sigmae(3)=curlADC(3)*curlADC(3)-Constant;
            sigmae(4)=curlADC(1)*curlADC(2);
            sigmae(5)=curlADC(1)*curlADC(3);
            sigmae(6)=curlADC(2)*curlADC(3);
            sigmae=(1/mue)*sigmae;
            
            ConstantD(1)=sigmae(1)*nm(1)+sigmae(4)*nm(2)+sigmae(5)*nm(3);
            ConstantD(2)=sigmae(4)*nm(1)+sigmae(2)*nm(2)+sigmae(6)*nm(3);
            ConstantD(3)=sigmae(5)*nm(1)+sigmae(6)*nm(2)+sigmae(3)*nm(3);


                        end
                    end
                end
            end
             
             
             
             R_U=R_U-Basis3Df*ConstantD*intfw(pp)*area;
             
      

              
            
 
          

 
end



