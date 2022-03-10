% Function to compute elemental boundary/interface integrals

function [K_UA,R_UR]= stiffLinearBounRes(K_UA,R_UR,xy,esize,esizeH1,intfxi,intfet,intfw,nipf,eltype,cond,...
    ielem,lec,lfc,flag,gorder,probdata,...
    mycoord,phh1f,phH1f,gphfx,gphfy,gphfz,ecphdfx,ecphdfy,ecphdfz,mue,probstatic,A,matc,subFlag,unknown,order,orderH1,mesh,X,DirichletFaceBasis)


ph=zeros(esize,3);
gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
gph=zeros(gesizet,3);
ph1=zeros(gesizet,1);
if flag==0
    gesizet=4;
end

globfa=mesh.face.globfa;
nelem=mesh.Nelements;
intma=mesh.intma;
coord=mesh.Coordinates;


% set up normals on the faces of the refeerence element
nmf(1,1)=-0.5;
nmf(1,2)=-(sqrt(3)/6);
nmf(1,3)=-(sqrt(6)/12);

nmf(2,1)=0.5;
nmf(2,2)=-(sqrt(3)/6);
nmf(2,3)=-(sqrt(6)/12);

nmf(3,1)=0;
nmf(3,2)=(sqrt(3)/3);
nmf(3,3)=-(sqrt(6)/12);

nmf(4,1)=0;
nmf(4,2)=0;
nmf(4,3)=(sqrt(3)/(2*sqrt(2)));

v(1,1)=-1;
v(1,2)=0;
v(1,3)=0;

v(2,1)=1;
v(2,2)=0;
v(2,3)=0;

v(3,1)=0;
v(3,2)=sqrt(3);
v(3,3)=0;

v(4,1)=0;
v(4,2)=sqrt(3)/3;
v(4,3)=2*(sqrt(2)/sqrt(3));

p=zeros(3,3);
contador=0;
for j=1:4
    if cond(ielem,j)==4 
        if matc(ielem)==1
              faceNumber=globfa(ielem,j);
              contador=contador+1;
              if contador>1
                  K_UA       = zeros(3*esizeH1,esize);
                  R_UR       = zeros(3*esizeH1,1);
              end
        %display('found bc face')
        % this is a boundary face
        % define normal's
        if j==1
            % face 1
            if eltype==1
                tau1(1)=v(3,1)-v(2,1);
                tau1(2)=v(3,2)-v(2,2);
                tau1(3)=v(3,3)-v(2,3);
                
                tau2(1)=v(4,1)-v(2,1);
                tau2(2)=v(4,2)-v(2,2);
                tau2(3)=v(4,3)-v(2,3);
            else
                tau1(1)=-v(3,1)+v(2,1);
                tau1(2)=-v(3,2)+v(2,2);
                tau1(3)=-v(3,3)+v(2,3);
                
                tau2(1)=v(4,1)-v(3,1);
                tau2(2)=v(4,2)-v(3,2);
                tau2(3)=v(4,3)-v(3,3);
            end
            % face 2
        elseif j==2
            tau1(1)=v(3,1)-v(1,1);
            tau1(2)=v(3,2)-v(1,2);
            tau1(3)=v(3,3)-v(1,3);
            
            tau2(1)=v(4,1)-v(1,1);
            tau2(2)=v(4,2)-v(1,2);
            tau2(3)=v(4,3)-v(1,3);
            
            % face 3
        elseif j==3
            tau1(1)=v(2,1)-v(1,1);
            tau1(2)=v(2,2)-v(1,2);
            tau1(3)=v(2,3)-v(1,3);
            
            tau2(1)=v(4,1)-v(1,1);
            tau2(2)=v(4,2)-v(1,2);
            tau2(3)=v(4,3)-v(1,3);
            
            % face 4
        else
            if eltype==1
                tau1(1)=v(2,1)-v(1,1);
                tau1(2)=v(2,2)-v(1,2);
                tau1(3)=v(2,3)-v(1,3);
                
                tau2(1)=v(3,1)-v(1,1);
                tau2(2)=v(3,2)-v(1,2);
                tau2(3)=v(3,3)-v(1,3);
                
            else
                tau1(1)=v(3,1)-v(1,1);
                tau1(2)=v(3,2)-v(1,2);
                tau1(3)=v(3,3)-v(1,3);
                
                tau2(1)=v(2,1)-v(1,1);
                tau2(2)=v(2,2)-v(1,2);
                tau2(3)=v(2,3)-v(1,3);
            end
        end
        
        % determine local face number and whether it is a boundary
        
        % nf contains local face number
        if j==1
            % set up integration point locations
            p(1,1:3)=v(2,1:3);
            p(2,1:3)=v(3,1:3);
            p(3,1:3)=v(4,1:3);
            
        elseif j==2
            p(1,1:3)=v(3,1:3);
            p(2,1:3)=v(1,1:3);
            p(3,1:3)=v(4,1:3);
            
        elseif j==3
            p(1,1:3)=v(1,1:3);
            p(2,1:3)=v(2,1:3);
            p(3,1:3)=v(4,1:3);
            
        else
            p(1,1:3)=v(1,1:3);
            p(2,1:3)=v(3,1:3);
            p(3,1:3)=v(2,1:3);
            
        end
        
        
        for pp=1:nipf
            S=zeros(3*esize,3);
             S2=zeros(3*esize,3);
            S_R=zeros(3,3);
            % compute integration point locations
            l(1)=1-intfxi(pp)-intfet(pp);
            l(2)=intfxi(pp);
            l(3)=intfet(pp);
            
            xi = l(1:3)*p(1:3,1);
            eta = l(1:3)*p(1:3,2);
            zeta = l(1:3)*p(1:3,3);
            gph(1:gesizet,1)=gphfx((j-1)*nipf+pp,1:gesizet)';
            gph(1:gesizet,2)=gphfy((j-1)*nipf+pp,1:gesizet)';
            gph(1:gesizet,3)=gphfz((j-1)*nipf+pp,1:gesizet)';
            ph1(1:gesizet,1)=phh1f((j-1)*nipf+pp,1:gesizet)';
            
            
            % evaluate covairant mapping
            [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord) ;
            
            %
            % evaluate basis
                    
                      Nhf=phH1f((j-1)*nipf+pp,1:esizeH1);
                     ph(1:esize,1:3)=(ecphdfx((j-1)*nipf+pp,1:esize)'*asxi(1:3))+...
                (ecphdfy((j-1)*nipf+pp,1:esize)'*aseta(1:3))+...
                (ecphdfz((j-1)*nipf+pp,1:esize)'*aszeta(1:3));
            
            % computx x,y,z
            [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
            
            nm(1)=(nmf(j,1)*axi(1))+(nmf(j,2)*aeta(1))+(nmf(j,3)*azeta(1));
            nm(2)=(nmf(j,1)*axi(2))+(nmf(j,2)*aeta(2))+(nmf(j,3)*azeta(2));
            nm(3)=(nmf(j,1)*axi(3))+(nmf(j,2)*aeta(3))+(nmf(j,3)*azeta(3));
            
            nd=sqrt((nm(1)^2)+(nm(2)^2)+(nm(3)^2));
            nm(1)=nm(1)/nd;
            nm(2)=nm(2)/nd;
            nm(3)=nm(3)/nd;
            
            
            % sign of computed normal is reversed from outward normal and so
            nm(1:3) = -nm(1:3);
            
            area=surdet(eltype,xi,eta,zeta,j,xy,lec,lfc,gorder+1);
            area=abs(area);
                
% 
% 

  for kk=1:nelem
                if matc(kk)==0
                    for ll=1:4
                        if globfa(kk,ll)==faceNumber
                             xy2 = coord(intma(kk,1:4),1:3);
                             flag2=0 ;
    gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
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
      
            
            [axi2,aeta2,azeta2,asxi2,aseta2,aszeta2,det2]=jacobian_pre(flag2,gesizet,gph,mycoord2) ;
            
                                 ph(1:esize,1:3)=(ecphdfx((ll-1)*nipf+pp,1:esize)'*asxi2(1:3))+...
                (ecphdfy((ll-1)*nipf+pp,1:esize)'*aseta2(1:3))+...
                (ecphdfz((ll-1)*nipf+pp,1:esize)'*aszeta2(1:3));
            
           curlADC=ph'*A;
           curlDeltaADC=ph'*A;
                 for ii=1:3
                    for b=1:esize
                        for jj=1:3
                            localrownumber2=(b-1)*3+ii;
                           S(localrownumber2,jj)=S(localrownumber2,jj)+(1/mue)*(curlADC(ii)*ph(b,jj)+ph(b,ii)*curlADC(jj)-ph(b,:)*curlADC*eq(ii,jj));
                           S2(b,ii,jj)=S2(localrownumber2,jj)+(1/mue)*(curlADC(ii)*ph(b,jj)*A(b)+ph(b,ii)*curlADC(jj)*A(b)-ph(b,:)*curlADC*eq(ii,jj)*A(b));
                        end
                    end
                 end
                   
                for ii=1:3
                        for jj=1:3
                           S_R(ii,jj)=S_R(ii,jj)+(1/mue)*(curlADC(ii)*curlDeltaADC(jj)+curlDeltaADC(ii)*curlADC(jj)-(curlDeltaADC'*curlADC)*eq(ii,jj));
                        end
                end
                        end
                    end
                end
            end
 

            
                for ii=1:3
                    for b=1:esize
                        localrownumber2=(b-1)*3+ii;
                                for a=1:esizeH1
                                 localrownumber=(a-1)*3+ii;
                                 localcolnumber=b;
                            for jj=1:3
                                 K_UA(localrownumber,localcolnumber)=K_UA(localrownumber,localcolnumber)-S(localrownumber2,jj)*(nm(jj))*Nhf(a)*area*intfw(pp);
                             
                            end
                                
              
                                end
                    end
                end
                
                
 

            
                for ii=1:3
                        
                                for a=1:esizeH1
                                 localrownumber=(a-1)*3+ii;
                            for jj=1:3
                                R_UR(localrownumber)=R_UR(localrownumber)-S_R(ii,jj)*(nm(jj))*Nhf(a)*area*intfw(pp);
                            end 
                                end
                end
          
        end
        end
    end
end
end



