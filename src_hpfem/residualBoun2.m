% Function to compute elemental boundary/interface integrals

function [R_U,AreaInt]= residualBoun2(R_U,xy,esize,esizeH1,intfxi,intfet,intfw,nipf,eltype,cond,...
    ielem,lec,lfc,flag,gorder,probdata,...
    mycoord,phh1f,phH1f,gphfx,gphfy,gphfz,ecphdfx,ecphdfy,ecphdfz,mue,probstatic,A,matc,AreaInt)


ph=zeros(esize,3);
gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
gph=zeros(gesizet,3);
ph1=zeros(gesizet,1);
if flag==0
    gesizet=4;
end



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


for j=1:4
    if cond(ielem,j)==4 
        if matc(ielem)==1
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
            sigmae=zeros(3,3);
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
            
            curlADC=ph'*A;
            curlADC=[1;0;0];

                for ii=1:3
                    for jj=1:3
                      sigmae(ii,jj)=(1/mue)*(curlADC(ii)*curlADC(jj)-0.5*(curlADC'*curlADC)*eq(ii,jj));
                    end
                end
              
                        for a=1:esizeH1
                            for ii=1:3
                                for jj=1:3
                          localrownumber=(a-1)*3+ii;
                          R_U(localrownumber)=R_U(localrownumber)-sigmae(ii,jj)*(nm(jj))*Nhf(a)*intfw(pp)*area;
                                end
                            end
                        end
              
%             

                AreaInt=AreaInt+intfw(pp)*area;
            
 
          
        end
        end
    end
    end
end



