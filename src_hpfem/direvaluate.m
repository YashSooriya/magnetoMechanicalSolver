function [ephdx1,ephdy1,ephdz1,ephdx2,ephdy2,ephdz2,ephdfx1,ephdfy1,ephdfz1,ephdfx2,ephdfy2,ephdfz2,...
          phh1f1,phh1f2,gphfx1,gphfy1,gphfz1,gphfx2,gphfy2,gphfz2,phH1e1,phH1e2,phH1f1,phH1f2,ecphdfx1,ecphdfy1,ecphdfz1,ecphdfx2,ecphdfy2,ecphdfz2]=direvaluate(nipe,xie,order,orderH1,...
          esizet,esizeH1,nipf,intfxi,intfet,gorder)

ephdx1 = zeros(6*nipe,esizet);
ephdy1 = zeros(6*nipe,esizet);
ephdz1 = zeros(6*nipe,esizet);
ephdx2 = zeros(6*nipe,esizet);
ephdy2 = zeros(6*nipe,esizet);
ephdz2 = zeros(6*nipe,esizet);

% compute basis functions on edges
for eltype=1:2
    for j=1:6
        if j==1
            if eltype==1
                xi1=1;
                et1=0;
                zt1=0;
                xi2=0;
                et2=sqrt(3);
                zt2=0;
            else
                xi2=1;
                et2=0;
                zt2=0;
                xi1=0;
                et1=sqrt(3);
                zt1=0;
            end
        elseif j==2
            xi1=-1;
            et1=0;
            zt1=0;
            xi2=0;
            et2=sqrt(3);
            zt2=0;
        elseif j==3
            xi1=-1;
            et1=0;
            zt1=0;
            xi2=1;
            et2=0;
            zt2=0;
        elseif j==4
            xi1=-1;
            et1=0;
            zt1=0;
            xi2=0;
            et2=sqrt(3)/3;
            zt2=2*(sqrt(2)/sqrt(3));
        elseif j==5
            xi1=1;
            et1=0;
            zt1=0;
            xi2=0;
            et2=sqrt(3)/3;
            zt2=2*(sqrt(2)/sqrt(3));
        else
            xi1=0;
            et1=sqrt(3);
            zt1=0;
            xi2=0;
            et2=sqrt(3)/3;
            zt2=2.*(sqrt(2)/sqrt(3));
        end
        
        for p=1:nipe
            
            if j==1
                xi=(0.5*(1-xie(p))*1)+(0.5*(1+xie(p))*0);
                eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*sqrt(3));
                zeta=0;
            elseif j==2
                xi=(0.5*(1.-xie(p))*(-1))+(0.5*(1+xie(p))*0);
                eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*sqrt(3));
                zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*0);
            elseif j==3
                xi=(0.5*(1.-xie(p))*(-1))+(0.5*(1+xie(p))*1);
                eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*0);
                zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*0);
            elseif j==4
                xi=(0.5*(1.-xie(p))*(-1))+(0.5*(1+xie(p))*0);
                eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*(sqrt(3.)/3));
                zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*2*(sqrt(2)/sqrt(3)));
            elseif j==5
                xi=(0.5*(1.-xie(p))*1)+(0.5*(1+xie(p))*0);
                eta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*(sqrt(3)/3));
                zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*2*(sqrt(2)/sqrt(3)));
            else
                xi=(0.5*(1.-xie(p))*0)+(0.5*(1+xie(p))*0);
                eta=(0.5*(1-xie(p))*sqrt(3))+(0.5*(1+xie(p))*(sqrt(3)/3));
                zeta=(0.5*(1-xie(p))*0)+(0.5*(1+xie(p))*2*(sqrt(2)/sqrt(3)));
            end
            
            % obtain basis
            % set up mapping functions
            axi(1)=1;
            axi(2)=0;
            axi(3)=0;
            
            aeta(1)=0;
            aeta(2)=1;
            aeta(3)=0;
            
            azeta(1)=0;
            azeta(2)=0;
            azeta(3)=1;
            
            asxi(1)=1;
            asxi(2)=0;
            asxi(3)=0;
    
            aseta(1)=0;
            aseta(2)=1;
            aseta(3)=0;
    
            aszeta(1)=0;
            aszeta(2)=0;
            aszeta(3)=1;
            
            ph= basis(xi,eta,zeta,axi,aeta,azeta,order,eltype,esizet);
            phH1=basish1(esizeH1,xi,eta,zeta,orderH1,eltype);
            
            if eltype==1
                % write(6,*)'evaluated basis',i
                ephdx1((j-1)*nipe+p,1:esizet)=ph(1:esizet,1)';
                ephdy1((j-1)*nipe+p,1:esizet)=ph(1:esizet,2)';
                ephdz1((j-1)*nipe+p,1:esizet)=ph(1:esizet,3)';
                phH1e1((j-1)*nipe+p,1:esizeH1)=phH1(1:esizeH1,1)';
                
            else
                ephdx2((j-1)*nipe+p,1:esizet)=ph(1:esizet,1)';
                ephdy2((j-1)*nipe+p,1:esizet)=ph(1:esizet,2)';
                ephdz2((j-1)*nipe+p,1:esizet)=ph(1:esizet,3)';
                phH1e2((j-1)*nipe+p,1:esizeH1)=phH1(1:esizeH1,1)';
                
            end
        end
    end
end


% compute basis functions on faces
ephdfx1=zeros(nipf*4,esizet);
ephdfy1=zeros(nipf*4,esizet);
ephdfz1=zeros(nipf*4,esizet);
ephdfx2=zeros(nipf*4,esizet);
ephdfy2=zeros(nipf*4,esizet);
ephdfz2=zeros(nipf*4,esizet);
ecphdfx1=zeros(nipf*4,esizet);
ecphdfy1=zeros(nipf*4,esizet);
ecphdfz1=zeros(nipf*4,esizet);
ecphdfx2=zeros(nipf*4,esizet);
ecphdfy2=zeros(nipf*4,esizet);
ecphdfz2=zeros(nipf*4,esizet);



% compute geometry basis functions on faces
gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
phh1f1=zeros(nipf*4,gesizet);
phh1f2=zeros(nipf*4,gesizet);
gphfx1=zeros(nipf*4,gesizet);
gphfy1=zeros(nipf*4,gesizet);
gphfz1=zeros(nipf*4,gesizet);
gphfx2=zeros(nipf*4,gesizet);
gphfy2=zeros(nipf*4,gesizet);
gphfz2=zeros(nipf*4,gesizet);



% vertices on reference element
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
for eltype=1:2
    for j=1:4
        if j==1
            % set up integration point locations
            pt(1,1:3)=v(2,1:3);
            pt(2,1:3)=v(3,1:3);
            pt(3,1:3)=v(4,1:3);
            
        elseif j==2
            pt(1,1:3)=v(3,1:3);
            pt(2,1:3)=v(1,1:3);
            pt(3,1:3)=v(4,1:3);
            
        elseif j==3
            pt(1,1:3)=v(1,1:3);
            pt(2,1:3)=v(2,1:3);
            pt(3,1:3)=v(4,1:3);
            
        else
            pt(1,1:3)=v(1,1:3);
            pt(2,1:3)=v(3,1:3);
            pt(3,1:3)=v(2,1:3);
        end
        
        
        for pp=1:nipf
            % write(6,*)intfxi(pp),intfet(pp),intfw(pp)
            % compute integration point locations
            l(1)=1-intfxi(pp)-intfet(pp);
            l(2)=intfxi(pp);
            l(3)=intfet(pp);
            
            xi = l(1:3)*pt(1:3,1);
            eta = l(1:3)*pt(1:3,2);
            zeta = l(1:3)*pt(1:3,3);
            
            % evaluate basis
            ph= basis(xi,eta,zeta,axi,aeta,azeta,order,eltype,esizet);
            
            if eltype==1
                % write(6,*)'evaluated basis',i
                ephdfx1((j-1)*nipf+pp,1:esizet)=ph(1:esizet,1)';
                ephdfy1((j-1)*nipf+pp,1:esizet)=ph(1:esizet,2)';
                ephdfz1((j-1)*nipf+pp,1:esizet)=ph(1:esizet,3)';
                
            else
                ephdfx2((j-1)*nipf+pp,1:esizet)=ph(1:esizet,1)';
                ephdfy2((j-1)*nipf+pp,1:esizet)=ph(1:esizet,2)';
                ephdfz2((j-1)*nipf+pp,1:esizet)=ph(1:esizet,3)';
                
            end
            
            % Evaluate curl basis
             ph= curlbasis(xi,eta,zeta,asxi,aseta,aszeta,order,eltype,esizet);
            
            if eltype==1
                % write(6,*)'evaluated basis',i
                ecphdfx1((j-1)*nipf+pp,1:esizet)=ph(1:esizet,1)';
                ecphdfy1((j-1)*nipf+pp,1:esizet)=ph(1:esizet,2)';
                ecphdfz1((j-1)*nipf+pp,1:esizet)=ph(1:esizet,3)';
                
            else
                ecphdfx2((j-1)*nipf+pp,1:esizet)=ph(1:esizet,1)';
                ecphdfy2((j-1)*nipf+pp,1:esizet)=ph(1:esizet,2)';
                ecphdfz2((j-1)*nipf+pp,1:esizet)=ph(1:esizet,3)';
                
            end
            
            % geometry information
            
            gph=gbasish1(xi,eta,zeta,axi,aeta,azeta,gorder+1,eltype,gesizet);
            
            if eltype==1
                gphfx1((j-1)*nipf+pp,1:gesizet)=gph(1:gesizet,1)';
                gphfy1((j-1)*nipf+pp,1:gesizet)=gph(1:gesizet,2)';
                gphfz1((j-1)*nipf+pp,1:gesizet)=gph(1:gesizet,3)';
            else
                gphfx2((j-1)*nipf+pp,1:gesizet)=gph(1:gesizet,1)';
                gphfy2((j-1)*nipf+pp,1:gesizet)=gph(1:gesizet,2)';
                gphfz2((j-1)*nipf+pp,1:gesizet)=gph(1:gesizet,3)';
            end
            
            phh1=basish1(gesizet,xi,eta,zeta,gorder+1,eltype);
            phH1=basish1(esizeH1,xi,eta,zeta,orderH1,eltype);
            if eltype==1
                phh1f1((j-1)*nipf+pp,1:gesizet)=phh1(1:gesizet,1)';
                phH1f1((j-1)*nipf+pp,1:esizeH1)=phH1(1:esizeH1,1)';
            else
                phh1f2((j-1)*nipf+pp,1:gesizet)=phh1(1:gesizet,1)';
                phH1f2((j-1)*nipf+pp,1:esizeH1)=phH1(1:esizeH1,1)';
            end
            
        end
    end
end

