function [hcurlerr,DispNorm] =ErrorDynamic(mesh,unknown,Basis,Quadrature,sol,probdata,probstatic,freq)

% Extract relevant data from mesh structure
nelem=mesh.Nelements;
eltype=mesh.eltype;
edgecof=mesh.edgecof;
facecof=mesh.facecof;
intma=mesh.intma;
glob=mesh.edge.glob;
globfa=mesh.face.globfa;
coord=mesh.Coordinates;
matc=mesh.matc;
mat=mesh.mat;


% Extract relevant data from ProblemData structure
esizet=probdata.esizet;
esizeH1=probdata.esizeH1;
gorder=probdata.jb.gorder;

% Extract relevant data from Basis structure
ephx1=Basis.ephx1;
ephy1=Basis.ephy1;
ephz1=Basis.ephz1;
ephx2=Basis.ephx2;
ephy2=Basis.ephy2;
ephz2=Basis.ephz2;
ecphx1=Basis.ecphx1;
ecphy1=Basis.ecphy1;
ecphz1=Basis.ecphz1;
ecphx2=Basis.ecphx2;
ecphy2=Basis.ecphy2;
ecphz2=Basis.ecphz2;
gphx1=Basis.gphx1;
gphy1=Basis.gphy1;
gphz1=Basis.gphz1;
gphx2=Basis.gphx2;
gphy2=Basis.gphy2;
gphz2=Basis.gphz2;
gphQuadx=Basis.gphQuadx;
gphQuady=Basis.gphQuady;
gphQuadz=Basis.gphQuadz;
intxi=Quadrature.intxi;
inteta=Quadrature.inteta;
intzeta=Quadrature.intzeta;
phh11=Basis.phh11;
phh12=Basis.phh12;
gradH1basis1x=Basis.gradH1basis1x;
gradH1basis1y=Basis.gradH1basis1y;
gradH1basis1z=Basis.gradH1basis1z;
gradH1basis2x=Basis.gradH1basis2x;
gradH1basis2y=Basis.gradH1basis2y;
gradH1basis2z=Basis.gradH1basis2z;

% Extract relevant data from Quadrature structure
intw=Quadrature.intw;
nip=Quadrature.nip;

probFlag=1;

lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);

error=0;
errord=0;
EnergyError=0;
ExactNormSum=0;
SNSerror=0;
exactNormSNSsum=0;
DispNorm=0;



% Define angular frequency and conductivity
omega=freq*pi*2;


for i=1:nelem
     if matc(i)==1
    xy = coord(intma(i,1:4),1:3);
    
    flag=0;
        gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
        mycoord=zeros(gesizet,3);
        mycoord(1:4,1:3)=xy;
       if gorder>1
        for j=1:6
            for pp=1:gorder
                for k=1:3
                    lec(j,pp,k)=edgecof(glob(i,j),((pp-1)*3)+k);
                    if(abs(edgecof(glob(i,j),((pp-1)*3)+k))>0.00000001)
                        flag=1;
                    end
                end
            end
        end
        
        for j=1:4
            for pp=1:gorder*(gorder-1)/2
                for k=1:3
                    lfc(j,pp,k)=facecof(globfa(i,j),((pp-1)*3)+k);
                end
            end
        end
        
        
        % edge functions
        for ii=1:6
            for p=1:gorder
                for j=1:3
                    mycoord(4+ii+6*(p-1),j)=lec(ii,p,j);
                end
            end
        end
        
        % face functions
        for iii=1:4
            for ii=1:(gorder-1)*gorder/2
                for j=1:3
                    mycoord(4+6*gorder+(iii-1)*gorder*(gorder-1)/2+ii,j)= lfc(iii,ii,j);
                end
            end
        end
        elseif gorder==1
            mycoord=mesh.mycoord(:,:,i);
            flag=1;
        end
        
        if flag==0
            gesizet=4;
        end
    
    %-------------------------------------------------------------------------

    lunkv=unknown.EM.unknowns(:,i);
    
    % extract local solution for this element
    lsol=zeros(esizet,1);
    lsol(lunkv>0)=sol(lunkv(lunkv>0));
    
    % choose correct set of stored basis functions
    if eltype(i)==1
        gphx=gphx1;
        gphy=gphy1;
        gphz=gphz1;
        phh1=phh11;
        ephx=ephx1;
        ephy=ephy1;
        ephz=ephz1;
        ecphx=ecphx1;
        ecphy=ecphy1;
        ecphz=ecphz1;
    else
        gphx=gphx2;
        gphy=gphy2;
        gphz=gphz2;
        phh1=phh12;
        ephx=ephx2;
        ephy=ephy2;
        ephz=ephz2;
        ecphx=ecphx2;
        ecphy=ecphy2;
        ecphz=ecphz2;
    end
    
    
    %-------------------------------------------------------------------------
    % evaluate covairant mapping (for linear geometry mapping is constant)
    if flag==0
        gph(1:gesizet,1)=gphx(1,1:gesizet)';
        gph(1:gesizet,2)=gphy(1,1:gesizet)';
        gph(1:gesizet,3)=gphz(1,1:gesizet)';
        
        [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
        
        
        for pp=1:nip
            
            ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
            % use stored functions
            ph(1:esizet,1:3)=(ephx(pp,1:esizet)'*axi(1:3))+...
                (ephy(pp,1:esizet)'*aeta(1:3))+...
                (ephz(pp,1:esizet)'*azeta(1:3)); 
            
            
            % compute the solution for this problem, at this integration point
            % in this element
            e=ph'*lsol; % This is the magnetic vector potential
            
            
            
            % use stored functions
            ph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                (ecphy(pp,1:esizet)'*aseta(1:3))+...
                (ecphz(pp,1:esizet)'*aszeta(1:3));
            
            % compute the solution for this problem, at this integration point
            % in this element
            curle=ph'*lsol;  % This is the magnetic flux density
            
            
            % computx x,y,z
            [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
            
            
            fun=probdata.es.exactfun;
            arg=probdata.es.exactfunarg;
            exact=fun(x,y,z,arg,probstatic,probFlag,omega);                        
            
            fun=probdata.es.exactcurlfun;
            arg=probdata.es.exactcurlfunarg;
            [curlexact,~,~]=fun(x,y,z,arg,probstatic,omega);
            
            for ii=1:3
                error=error+intw(pp)*det*(abs(e(ii)-exact(ii))^2+...
                    abs(curle(ii)-curlexact(ii))^2);
%                                 error=error+intw(pp)*det*(abs(e(ii)-exact(ii))^2);
                errord=errord+intw(pp)*det*(abs(exact(ii))^2+abs(curlexact(ii))^2);


            end
            
      
            
            
        end
        
    else
        
        for pp=1:nip
                if gorder==1
                gphQuad(1:gesizet,1)=gphQuadx(pp,1:gesizet)';
                gphQuad(1:gesizet,2)=gphQuady(pp,1:gesizet)';
                gphQuad(1:gesizet,3)=gphQuadz(pp,1:gesizet)';
                
                
                % evaluate covairant mapping
                [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gphQuad,mycoord);
                [x,y,z]= getxyzq(mycoord,intxi(pp),inteta(pp),intzeta(pp));
                else
                gph(1:gesizet,1)=gphx(pp,1:gesizet)';
                gph(1:gesizet,2)=gphy(pp,1:gesizet)';
                gph(1:gesizet,3)=gphz(pp,1:gesizet)';
                ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
                
                % evaluate covairant mapping
                [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
                [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);  
                end
            
            % use stored functions            
            ph(1:esizet,1:3)=(ephx(pp,1:esizet)'*axi(1:3))+...
                (ephy(pp,1:esizet)'*aeta(1:3))+...
                (ephz(pp,1:esizet)'*azeta(1:3));
                  
            % compute the solution for this problem, at this integration point
            % in this element          
            e=ph'*lsol;
            
            
            %             % use stored functions
            ph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                (ecphy(pp,1:esizet)'*aseta(1:3))+...
                (ecphz(pp,1:esizet)'*aszeta(1:3));
            
            % compute the solution for this problem, at this integration point
            % in this element
            curle=ph'*lsol;
            
            
            fun=probdata.es.exactfun;
            arg=probdata.es.exactfunarg;
            exact=fun(x,y,z,arg,probstatic,probFlag,omega);
            
            
            
            fun=probdata.es.exactcurlfun;
            arg=probdata.es.exactcurlfunarg;
            [curlexact,~,~]=fun(x,y,z,arg,probstatic,omega);
            
            for ii=1:3
                error=error+intw(pp)*det*(abs(e(ii)-exact(ii))^2+...
                    abs(curle(ii)-curlexact(ii))^2);

                errord=errord+intw(pp)*det*(abs(exact(ii))^2+abs(curlexact(ii))^2);
            end
            
            
            
            
        end
        
    end
     end
    % end of loop over elements
end
hcurlerr = sqrt(error/errord);
display(['The hcurl error is = ',num2str(hcurlerr)])

%==================================================================================================
% Calculate error for the mechanical problem
%==================================================================================================

% Extract relevant data from mesh structure
nelem=mesh.mech.Nelements;
eltype=mesh.mech.eltype;
edgecof=mesh.mech.edgecof;
facecof=mesh.mech.facecof;
intma=mesh.mech.intma;
glob=mesh.mech.glob;
globfa=mesh.mech.globfa;
coord=mesh.mech.Coordinates;
mapL2G_e=mesh.mech.mapL2G_e;

% Extract relevant data from Basis structure
gphx1=Basis.gphx1;
gphy1=Basis.gphy1;
gphz1=Basis.gphz1;
gphx2=Basis.gphx2;
gphy2=Basis.gphy2;
gphz2=Basis.gphz2;
phh11=Basis.phh11;
phh12=Basis.phh12;
H1basis1=Basis.H1basis1;
H1basis2=Basis.H1basis2;

% Extract relevant data from Quadrature structure
intw=Quadrature.intw;
nip=Quadrature.nip;

% Define problem flag (2 for mechanics)
probFlag=2;

% Initialise variables
lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);

error=0;
errord=0;


%--------------------------------------------------------------------------
% Loop over elements
%--------------------------------------------------------------------------
for i=1:nelem
    % Local to global mapping
    aaa=mapL2G_e(i);
    material=mat(i);
    % Extract Lame parameters from structure
    lambda=probdata.matr.lambda(material);
    G=probdata.matr.G(material);
    
    % Define coordinates of vertices
    xy = coord(intma(i,1:4),1:3);
    
    flag=0;
        gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
        mycoord=zeros(gesizet,3);
        mycoord(1:4,1:3)=xy;
       if gorder>1
        for j=1:6
            for pp=1:gorder
                for k=1:3
                    lec(j,pp,k)=edgecof(glob(i,j),((pp-1)*3)+k);
                    if(abs(edgecof(glob(i,j),((pp-1)*3)+k))>0.00000001)
                        flag=1;
                    end
                end
            end
        end
        
        for j=1:4
            for pp=1:gorder*(gorder-1)/2
                for k=1:3
                    lfc(j,pp,k)=facecof(globfa(i,j),((pp-1)*3)+k);
                end
            end
        end

        
        % edge functions
        for ii=1:6
            for p=1:gorder
                for j=1:3
                    mycoord(4+ii+6*(p-1),j)=lec(ii,p,j);
                end
            end
        end
        
        % face functions
        for iii=1:4
            for ii=1:(gorder-1)*gorder/2
                for j=1:3
                    mycoord(4+6*gorder+(iii-1)*gorder*(gorder-1)/2+ii,j)= lfc(iii,ii,j);
                end
            end
        end
        elseif gorder==1
            mycoord=mesh.mycoord(:,:,aaa);
            flag=1;
        end
        
        if flag==0
            gesizet=4;
        end
    %-------------------------------------------------------------------------
    % work out numbering of basis functions
    lunkv=unknown.Mech.unknownsD(:,aaa);
    
    
    % extract local solution for this element
    lsol=zeros(3*esizeH1,1);
    lsol(lunkv>0)=sol(lunkv(lunkv>0));
    
    % choose correct set of stored basis functions
    if eltype(i)==1
        gphx=gphx1;
        gphy=gphy1;
        gphz=gphz1;
        phh1=phh11;
        H1basis=H1basis1;
        gradH1basisx=gradH1basis1x;
        gradH1basisy=gradH1basis1y;
        gradH1basisz=gradH1basis1z;
    else
        gphx=gphx2;
        gphy=gphy2;
        gphz=gphz2;
        phh1=phh12;
        H1basis=H1basis2;
        gradH1basisx=gradH1basis2x;
        gradH1basisy=gradH1basis2y;
        gradH1basisz=gradH1basis2z;
    end
    
    
    %-------------------------------------------------------------------------
    % evaluate covairant mapping (for linear geometry mapping is constant)
    if flag==0
        gph(1:gesizet,1)=gphx(1,1:gesizet)';
        gph(1:gesizet,2)=gphy(1,1:gesizet)';
        gph(1:gesizet,3)=gphz(1,1:gesizet)';
        
        [axi,aeta,azeta,~,~,~,det]=jacobian_pre(flag,gesizet,gph,mycoord);
        
        
        for pp=1:nip
            
            ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
            
            % use stored functions
            
            H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
            H1bas3D=zeros(3,3*esizeH1);
            gradH1bas3D_x=zeros(3,3*esizeH1);
            gradH1bas3D_y=zeros(3,3*esizeH1);
            gradH1bas3D_z=zeros(3,3*esizeH1);
            
            H1bas3D(1,1:3:end)=H1bas;
            H1bas3D(2,2:3:end)=H1bas;
            H1bas3D(3,3:3:end)=H1bas;
            
            % Construct matrices to calculate energy error
            
            gradH1bas(1:esizeH1,1:3)=(gradH1basisx(pp,1:esizeH1)'*axi(1:3))+...
                (gradH1basisy(pp,1:esizeH1)'*aeta(1:3))+...
                (gradH1basisz(pp,1:esizeH1)'*azeta(1:3));
            
            gradH1bas3D_x(1,1:3:end)=gradH1bas(:,1);
            gradH1bas3D_x(2,1:3:end)=gradH1bas(:,2);
            gradH1bas3D_x(3,1:3:end)=gradH1bas(:,3);
            gradH1bas3D_y(1,2:3:end)=gradH1bas(:,1);
            gradH1bas3D_y(2,2:3:end)=gradH1bas(:,2);
            gradH1bas3D_y(3,2:3:end)=gradH1bas(:,3);
            gradH1bas3D_z(1,3:3:end)=gradH1bas(:,1);
            gradH1bas3D_z(2,3:3:end)=gradH1bas(:,2);
            gradH1bas3D_z(3,3:3:end)=gradH1bas(:,3);
            
            
            a_1=gradH1bas3D_x*lsol;
            a_2=gradH1bas3D_y*lsol;
            a_3=gradH1bas3D_z*lsol;
            
            grad_hp=[a_1';a_2';a_3'];
            %grad_hp=abs(grad_hp);
            div_hp=trace(grad_hp);
            strain_hp=0.5*(grad_hp+grad_hp');
            
            
            % compute the solution for this problem, at this integration point
            % in this element
            u_D=zeros(3,1);
            u_D=H1bas3D*lsol;
            
            
            
            
            
            % computx x,y,z
            [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
            
            
            fun=probdata.es.exactfun;
            arg=probdata.es.exactfunarg;
            exact=fun(x,y,z,arg,probstatic,probFlag);
            
            func=probdata.es.GradientExact;
            grad_exact=func(x,y,z,arg);
            div_exact=trace(grad_exact);
            strain_exact=0.5*(grad_exact+grad_exact');
            
            strain_difference=strain_exact-strain_hp;
            strain_norm=(strain_difference(1,1)^2+strain_difference(1,2)^2+strain_difference(1,3)^2+strain_difference(2,1)^2+strain_difference(2,2)^2+strain_difference(2,3)^2+...
                +strain_difference(3,1)^2+strain_difference(3,2)^2+strain_difference(3,3)^2);
            div_norm=abs(div_exact-div_hp)^2;
            errorNorm=(2*G*strain_norm+lambda*div_norm);
            
            strain_exact_norm=(strain_exact(1,1)^2+strain_exact(1,2)^2+strain_exact(1,3)^2+strain_exact(2,1)^2+strain_exact(2,2)^2+strain_exact(2,3)^2+...
                +strain_exact(3,1)^2+strain_exact(3,2)^2+strain_exact(3,3)^2);
            exactNorm=(2*G*strain_exact_norm+lambda*div_exact^2);
            
            
            
            EnergyError=EnergyError+intw(pp)*det*errorNorm;
            ExactNormSum=ExactNormSum+intw(pp)*det*exactNorm;
            
            % Now construct the stress tensor to calculate SNS error
            
            Idmat=zeros(3,3);
            Idmat(1,1)=1;
            Idmat(2,2)=1;
            Idmat(3,3)=1;
            
            % Construct the exact stress tensor
            
            StressTensor_exact=lambda*div_exact*Idmat+2*G*strain_exact;
            StressTensor_hp=lambda*div_hp*Idmat+2*G*strain_hp;
            
            errorNormSNS=(StressTensor_exact(1,1)+StressTensor_exact(2,2)+StressTensor_exact(3,3)-(StressTensor_hp(1,1)+StressTensor_hp(2,2)+StressTensor_hp(3,3)))^2;
            exactNormSNS=(StressTensor_exact(1,1)+StressTensor_exact(2,2)+StressTensor_exact(3,3))^2;
            
            SNSerror=SNSerror+intw(pp)*det*errorNormSNS;
            exactNormSNSsum=exactNormSNSsum+exactNormSNS*intw(pp)*det;
            
            
            
            
            
            
            for ii=1:3
                error=error+intw(pp)*det*(abs(u_D(ii)-exact(ii)))^2;
                
                
                errord=errord+intw(pp)*det*(abs(exact(ii)))^2;
                DispNorm=DispNorm+intw(pp)*det*norm(u_D(ii))*norm(u_D(ii));
                
            end
            
            
        end
        
    else
        
        for pp=1:nip
                if gorder==1
                gphQuad(1:gesizet,1)=gphQuadx(pp,1:gesizet)';
                gphQuad(1:gesizet,2)=gphQuady(pp,1:gesizet)';
                gphQuad(1:gesizet,3)=gphQuadz(pp,1:gesizet)';
                
                
                % evaluate covairant mapping
                [axi,aeta,azeta,~,~,~,det]=jacobian_pre(flag,gesizet,gphQuad,mycoord);
                [x,y,z]= getxyzq(mycoord,intxi(pp),inteta(pp),intzeta(pp));
                else
                gph(1:gesizet,1)=gphx(pp,1:gesizet)';
                gph(1:gesizet,2)=gphy(pp,1:gesizet)';
                gph(1:gesizet,3)=gphz(pp,1:gesizet)';
                ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
                
                % evaluate covairant mapping
                [axi,aeta,azeta,~,~,~,det]=jacobian_pre(flag,gesizet,gph,mycoord);
                [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);  
                end
            
            % use stored functions
            
            H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
            H1bas3D=zeros(3,3*esizeH1);
            gradH1bas3D_x=zeros(3,3*esizeH1);
            gradH1bas3D_y=zeros(3,3*esizeH1);
            gradH1bas3D_z=zeros(3,3*esizeH1);
            
            H1bas3D(1,1:3:end)=H1bas;
            H1bas3D(2,2:3:end)=H1bas;
            H1bas3D(3,3:3:end)=H1bas;
            
            % Construct matrices to calculate energy error
            
            gradH1bas(1:esizeH1,1:3)=(gradH1basisx(pp,1:esizeH1)'*axi(1:3))+...
                (gradH1basisy(pp,1:esizeH1)'*aeta(1:3))+...
                (gradH1basisz(pp,1:esizeH1)'*azeta(1:3));
            
            gradH1bas3D_x(1,1:3:end)=gradH1bas(:,1);
            gradH1bas3D_x(2,1:3:end)=gradH1bas(:,2);
            gradH1bas3D_x(3,1:3:end)=gradH1bas(:,3);
            gradH1bas3D_y(1,2:3:end)=gradH1bas(:,1);
            gradH1bas3D_y(2,2:3:end)=gradH1bas(:,2);
            gradH1bas3D_y(3,2:3:end)=gradH1bas(:,3);
            gradH1bas3D_z(1,3:3:end)=gradH1bas(:,1);
            gradH1bas3D_z(2,3:3:end)=gradH1bas(:,2);
            gradH1bas3D_z(3,3:3:end)=gradH1bas(:,3);
            
            
            a_1=gradH1bas3D_x*lsol;
            a_2=gradH1bas3D_y*lsol;
            a_3=gradH1bas3D_z*lsol;
            
            grad_hp=[a_1';a_2';a_3'];
            div_hp=trace(grad_hp);
            strain_hp=0.5*(grad_hp+grad_hp');
            
            
            % compute the solution for this problem, at this integration point
            % in this element
            u_D=H1bas3D*lsol;
            

            fun=probdata.es.exactfun;
            arg=probdata.es.exactfunarg;
            exact=fun(x,y,z,arg,probstatic,probFlag);
            func=probdata.es.GradientExact;
            grad_exact=func(x,y,z,arg);
            
            div_exact=trace(grad_exact);
            strain_exact=0.5*(grad_exact+grad_exact');
            
            strain_difference=strain_exact-strain_hp;
            strain_norm=(strain_difference(1,1)^2+strain_difference(1,2)^2+strain_difference(1,3)^2+strain_difference(2,1)^2+strain_difference(2,2)^2+strain_difference(2,3)^2+...
                +strain_difference(3,1)^2+strain_difference(3,2)^2+strain_difference(3,3)^2);
            div_norm=abs(div_exact-div_hp)^2;
            errorNorm=(2*G*strain_norm+lambda*div_norm);
            
            strain_exact_norm=(strain_exact(1,1)^2+strain_exact(1,2)^2+strain_exact(1,3)^2+strain_exact(2,1)^2+strain_exact(2,2)^2+strain_exact(2,3)^2+...
                +strain_exact(3,1)^2+strain_exact(3,2)^2+strain_exact(3,3)^2);
            exactNorm=(2*G*strain_exact_norm+lambda*div_exact^2);
            
            
            EnergyError=EnergyError+intw(pp)*det*errorNorm;
            ExactNormSum=ExactNormSum+intw(pp)*det*exactNorm;
            
            
            
            % Now construct the stress tensor to calculate SNS error
            
            Idmat=zeros(3,3);
            Idmat(1,1)=1;
            Idmat(2,2)=1;
            Idmat(3,3)=1;
            
            % Construct the exact stress tensor
            
            StressTensor_exact=lambda*div_exact*Idmat+2*G*strain_exact;
            StressTensor_hp=lambda*div_hp*Idmat+2*G*strain_hp;
            
            errorNormSNS=(StressTensor_exact(1,1)+StressTensor_exact(2,2)+StressTensor_exact(3,3)-(StressTensor_hp(1,1)+StressTensor_hp(2,2)+StressTensor_hp(3,3)))^2;
            exactNormSNS=(StressTensor_exact(1,1)+StressTensor_exact(2,2)+StressTensor_exact(3,3))^2;
            
            SNSerror=SNSerror+intw(pp)*det*errorNormSNS;
            exactNormSNSsum=exactNormSNSsum+exactNormSNS*intw(pp)*det;
            
            
            
            
            for ii=1:3
                error=error+intw(pp)*det*(abs(u_D(ii)-exact(ii)))^2;
                DispNorm=DispNorm+intw(pp)*det*norm(u_D(ii))*norm(u_D(ii));
                
                
                errord=errord+intw(pp)*det*(abs(exact(ii)))^2;
                
                
            end
            
            
        end
        
    end
    
    % end of loop over elements
end
L2err_mech = sqrt(error/errord);
Energy_error_mech=abs(sqrt(EnergyError/ExactNormSum));
SNS_error_mech=abs(sqrt(SNSerror/exactNormSNSsum));
DispNorm=sqrt(DispNorm);
display(['The Displacement norm is = ',num2str(DispNorm)])
display(['The L2 error is = ',num2str(L2err_mech)])
display(['The energy error is = ',num2str(Energy_error_mech)])
display(['The SNS error is = ',num2str(SNS_error_mech)])