function [Basis]= evaluate(Quadrature,ProblemData)

% Extract relevant data from quadrature structure
intxi=Quadrature.intxi;
inteta=Quadrature.inteta;
intzeta=Quadrature.intzeta;
nip=Quadrature.nip;

% Extract relevant data from FEspacesInfo structure
order=ProblemData.order;
orderH1=ProblemData.orderH1;
esizet=ProblemData.esizet;
esizeH1=ProblemData.esizeH1;
gorder=ProblemData.jb.gorder;


ephx1 = zeros(nip,esizet);
ephy1 = zeros(nip,esizet);
ephz1 = zeros(nip,esizet);

ephx2 = zeros(nip,esizet);
ephy2 = zeros(nip,esizet);
ephz2 = zeros(nip,esizet);

ecphx1 = zeros(nip,esizet);
ecphy1 = zeros(nip,esizet);
ecphz1 = zeros(nip,esizet);

ecphx2 = zeros(nip,esizet);
ecphy2 = zeros(nip,esizet);
ecphz2 = zeros(nip,esizet);

gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;

gphx1=zeros(nip,gesizet);
gphy1=zeros(nip,gesizet);
gphz1=zeros(nip,gesizet);

gphx2=zeros(nip,gesizet);
gphy2=zeros(nip,gesizet);
gphz2=zeros(nip,gesizet);

phh11=zeros(nip,gesizet);
phh12=zeros(nip,gesizet);

gradH1basis1x=zeros(nip,esizeH1);
gradH1basis1y=zeros(nip,esizeH1);
gradH1basis1z=zeros(nip,esizeH1);

gradH1basis2x=zeros(nip,esizeH1);
gradH1basis2y=zeros(nip,esizeH1);
gradH1basis2z=zeros(nip,esizeH1);


for i=1:nip
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
    % write(6,*)intxi(i),inteta(i),intzeta(i),i,nip
    ph= basis(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,order,1,esizet);
    
    % write(6,*)'evaluated basis',i
    ephx1(i,1:esizet)=ph(1:esizet,1)';
    ephy1(i,1:esizet)=ph(1:esizet,2)';
    ephz1(i,1:esizet)=ph(1:esizet,3)';
    
    ph= basis(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,order,2,esizet);
    
    % write(6,*)'evaluated basis typ 2',i
    ephx2(i,1:esizet)=ph(1:esizet,1)';
    ephy2(i,1:esizet)=ph(1:esizet,2)';
    ephz2(i,1:esizet)=ph(1:esizet,3)';
    
    ph= curlbasis(intxi(i),inteta(i),intzeta(i),asxi,aseta,aszeta,order,1,esizet);
    ecphx1(i,1:esizet)=ph(1:esizet,1)';
    ecphy1(i,1:esizet)=ph(1:esizet,2)';
    ecphz1(i,1:esizet)=ph(1:esizet,3)';
    
    ph= curlbasis(intxi(i),inteta(i),intzeta(i),asxi,aseta,aszeta,order,2,esizet);
    ecphx2(i,1:esizet)=ph(1:esizet,1)';
    ecphy2(i,1:esizet)=ph(1:esizet,2)';
    ecphz2(i,1:esizet)=ph(1:esizet,3)';
    
    
    gph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,gorder+1,1,gesizet);
    
    gphx1(i,1:gesizet)=gph(1:gesizet,1)';
    gphy1(i,1:gesizet)=gph(1:gesizet,2)';
    gphz1(i,1:gesizet)=gph(1:gesizet,3)';
    
    gph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,gorder+1,2,gesizet);
    gphx2(i,1:gesizet)=gph(1:gesizet,1)';
    gphy2(i,1:gesizet)=gph(1:gesizet,2)';
    gphz2(i,1:gesizet)=gph(1:gesizet,3)';
    
    
    phh1=basish1(gesizet,intxi(i),inteta(i),intzeta(i),gorder+1,1);
    phh11(i,1:gesizet)=phh1(1:gesizet,1)';
    
    phh1=basish1(gesizet,intxi(i),inteta(i),intzeta(i),gorder+1,2);
    phh12(i,1:gesizet)=phh1(1:gesizet,1)';
    
    
    %Evaluate H1 basis functions for polynomial order p
    th1=basish1(esizeH1,intxi(i),inteta(i),intzeta(i),orderH1,1);
    H1basis1(i,1:esizeH1)=th1(1:esizeH1,1)';
    
    th1=basish1(esizeH1,intxi(i),inteta(i),intzeta(i),orderH1,2);
    H1basis2(i,1:esizeH1)=th1(1:esizeH1,1)';
    
    % Evaluate gradient of H1 basis functions for polynomial order p
    
    tgph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,orderH1,1,esizeH1);
    
    gradH1basis1x(i,1:esizeH1)=tgph(1:esizeH1,1)';
    gradH1basis1y(i,1:esizeH1)=tgph(1:esizeH1,2)';
    gradH1basis1z(i,1:esizeH1)=tgph(1:esizeH1,3)';
    
    tgph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,orderH1,2,esizeH1);
    gradH1basis2x(i,1:esizeH1)=tgph(1:esizeH1,1)';
    gradH1basis2y(i,1:esizeH1)=tgph(1:esizeH1,2)';
    gradH1basis2z(i,1:esizeH1)=tgph(1:esizeH1,3)';
    
    
end


% Store basis functions to Basis structure
Basis.ephx1=ephx1;
Basis.ephy1=ephy1;
Basis.ephz1=ephz1;
Basis.ephx2=ephx2;
Basis.ephy2=ephy2;
Basis.ephz2=ephz2;
Basis.ecphx1=ecphx1;
Basis.ecphy1=ecphy1;
Basis.ecphz1=ecphz1;
Basis.ecphx2=ecphx2;
Basis.ecphy2=ecphy2;
Basis.ecphz2=ecphz2;
Basis.gphx1=gphx1;
Basis.gphy1=gphy1;
Basis.gphz1=gphz1;
Basis.gphx2=gphx2;
Basis.gphy2=gphy2;
Basis.gphz2=gphz2;
Basis.phh11=phh11;
Basis.phh12=phh12;
Basis.H1basis1=H1basis1;
Basis.H1basis2=H1basis2;
Basis.gradH1basis1x=gradH1basis1x;
Basis.gradH1basis1y=gradH1basis1y;
Basis.gradH1basis1z=gradH1basis1z;
Basis.gradH1basis2x=gradH1basis2x;
Basis.gradH1basis2y=gradH1basis2y;
Basis.gradH1basis2z=gradH1basis2z;


%========================================================================
% Compute basis functions on edges and faces to apply boundary conditions
%========================================================================

% Extract relevant data from Quadrature structure
intfxi=Quadrature.intfxi;
intfet=Quadrature.intfet;
intfw=Quadrature.intfw;
nipf=Quadrature.nipf;

% Extract relevant data from ProblemData structure
order=ProblemData.order;
orderH1=ProblemData.orderH1;
esizet=ProblemData.esizet;
esizeH1=ProblemData.esizeH1;
gorder=ProblemData.jb.gorder;


% find integration points along edge
% int_-1^+1 dx
kind=1;
nipe=floor((2*(order+1)+1)/2)+1;
if nipe>20
    error(message('increase nipe in pec.m'));
end
kpts=0;
endpts(1)=-1;
endpts(2)=1;
[bb,xie,w]=gaussq1(kind,nipe,0,0,kpts,endpts);

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

Basis.Edge.ephdx1=ephdx1;
Basis.Edge.ephdy1=ephdy1;
Basis.Edge.ephdz1=ephdz1;
Basis.Edge.ephdx2=ephdx2;
Basis.Edge.ephdy2=ephdy2;
Basis.Edge.ephdz2=ephdz2;
Basis.Edge.phH1e1=phH1e1;
Basis.Edge.phH1e2=phH1e2;

Basis.Face.ephdfx1=ephdfx1;
Basis.Face.ephdfy1=ephdfy1;
Basis.Face.ephdfz1=ephdfz1;
Basis.Face.ephdfx2=ephdfx2;
Basis.Face.ephdfy2=ephdfy2;
Basis.Face.ephdfz2=ephdfz2;

Basis.Face.phh1f1=phh1f1;
Basis.Face.phh1f2=phh1f2;
Basis.Face.gphfx1=gphfx1;
Basis.Face.gphfy1=gphfy1;
Basis.Face.gphfz1=gphfz1;
Basis.Face.gphfx2=gphfx2;
Basis.Face.gphfy2=gphfy2;
Basis.Face.gphfz2=gphfz2;
Basis.Face.phH1f1=phH1f1;
Basis.Face.phH1f2=phH1f2;

Basis.Face.ecphdfx1=ecphdfx1;
Basis.Face.ecphdfy1=ecphdfy1;
Basis.Face.ecphdfz1=ecphdfz1;
Basis.Face.ecphdfx2=ecphdfx2;
Basis.Face.ecphdfy2=ecphdfy2;
Basis.Face.ecphdfz2=ecphdfz2;




