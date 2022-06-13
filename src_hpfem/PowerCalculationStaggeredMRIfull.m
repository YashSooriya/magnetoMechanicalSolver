function [OutPower4Kfull,OutPower77Kfull,OutPowerOVCfull] = PowerCalculationStaggeredMRIfull(Mesh,Unknown,Basis,Quadrature,...
    ProblemData,UnknownStatic,solStatic,condOVC,cond77K,cond4K,Dynamic,freqOut)


% Extract relevant data from mesh structure
nelem=Mesh.Nelements;
eltype=Mesh.eltype;
edgecof=Mesh.edgecof;
facecof=Mesh.facecof;
intma=Mesh.intma;
glob=Mesh.edge.glob;
globfa=Mesh.face.globfa;
coord=Mesh.Coordinates;
mat=Mesh.mat;
subFlag=Mesh.subFlag;

% Extract relevant data from ProblemData structure
esizet=ProblemData.esizet;
gorder=ProblemData.jb.gorder;
esizeH1=ProblemData.esizeH1;
matOVC=ProblemData.matOVC;
mat77K=ProblemData.mat77K;
mat4K=ProblemData.mat4K;

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
H1basis1=Basis.H1basis1;
H1basis2=Basis.H1basis2;

% Extract relevant data from Quadrature structure
intw=Quadrature.intw;
nip=Quadrature.nip;


lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);

OutPower4Kfull=zeros(length(freqOut),1);
OutPower77Kfull=zeros(length(freqOut),1);
OutPowerOVCfull=zeros(length(freqOut),1);
OutPower4Kfulltest=zeros(length(freqOut),1);

% Define angular frequency and conductivity
omegafull = zeros(length(freqOut),1);
for freqs=1:length(freqOut)
    omegafull(freqs)=freqOut(freqs)*pi*2;
end
sigma=ProblemData.matr.sigma;


for i=1:nelem
    if subFlag(i)==2
        xy = coord(intma(i,1:4),1:3);
        material=mat(i);
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
            mycoord=Mesh.mycoord(:,:,i);
            flag=1;
        end
        
        if flag==0
            gesizet=4;
        end
        
        %-----------------------------------------------------------------
        % Work out numbering of basis functions
        %-----------------------------------------------------------------
        lunkv=Unknown.EM.unknowns(:,i);
        lunkv2=Unknown.Mech.unknownsD(:,i);
        lunkvStatic=UnknownStatic.EM.unknowns(:,i);
        
        % Extract local solution for this element
        ldyn = zeros(esizet, length(freqOut));
        ldynmech = zeros(3*esizeH1, length(freqOut));
        
        for freqs = 1:length(freqOut)
            ldyn(:,freqs) = Dynamic(lunkv,freqs);
            ldynmech(:,freqs) = Dynamic(lunkv2,freqs);
        end

        lsolStatic=zeros(esizet,1);
        lsolStatic(lunkvStatic>0)=solStatic(lunkvStatic(lunkvStatic>0),1);
        
        % choose correct set of stored basis functions
        if eltype(i)==1
            gphx=gphx1;
            gphy=gphy1;
            gphz=gphz1;
            ephx=ephx1;
            ephy=ephy1;
            ephz=ephz1;
            ecphx=ecphx1;
            ecphy=ecphy1;
            ecphz=ecphz1;
            H1basis=H1basis1;
        else
            gphx=gphx2;
            gphy=gphy2;
            gphz=gphz2;
            ephx=ephx2;
            ephy=ephy2;
            ephz=ephz2;
            ecphx=ecphx2;
            ecphy=ecphy2;
            ecphz=ecphz2;
            H1basis=H1basis2;
        end
        
        
        %-------------------------------------------------------------------------
        % evaluate covairant mapping (for linear geometry mapping is constant)
        if flag==0
            gph(1:gesizet,1)=gphx(1,1:gesizet)';
            gph(1:gesizet,2)=gphy(1,1:gesizet)';
            gph(1:gesizet,3)=gphz(1,1:gesizet)';
            
            [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
            
            
            for pp=1:nip
                
                % use stored functions
                ph(1:esizet,1:3)=(ephx(pp,1:esizet)'*axi(1:3))+...
                    (ephy(pp,1:esizet)'*aeta(1:3))+...
                    (ephz(pp,1:esizet)'*azeta(1:3));
                
                  cph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                    (ecphy(pp,1:esizet)'*aseta(1:3))+...
                    (ecphz(pp,1:esizet)'*aszeta(1:3));
                
                
                % Compute the magnetic vector Potential
                efull = ph'*ldyn; % This is the magnetic vector potential
                
                % Compute the static magnetic flux density
                curlADC=cph'*lsolStatic;
            
                H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
                H1bas3D=zeros(3,3*esizeH1);

            
                H1bas3D(1,1:3:end)=H1bas;
                H1bas3D(2,2:3:end)=H1bas;
                H1bas3D(3,3:3:end)=H1bas;

            
                % compute the solution for this problem, at this integration point
                % in this element 
                u_Dfull = H1bas3D*ldynmech;
                

                
                % Compute the electric field                
                curlACDrep = repmat(curlADC, 1, length(freqOut));
                crossfull = cross(curlACDrep,u_Dfull);
                ElectricFieldfull = (1i*(omegafull.')).*((crossfull)-(efull));
                EF_full = abs(diag(ElectricFieldfull'*ElectricFieldfull));
                
                % The output power is:
       
%                 for freqs = 1:length(freqOut)
%                     EF = ElectricFieldfull(:,freqs);
%                     if any(material==mat4K)
%                         OutPower4Kfull(freqs)=OutPower4Kfull(freqs)+0.5*cond4K*sigma(material)*abs((EF'*EF))*intw(pp)*det;
%                     elseif any(material==mat77K)
%                         OutPower77Kfull(freqs)=OutPower77Kfull(freqs)+0.5*cond77K*sigma(material)*abs((EF'*EF))*intw(pp)*det;
%                     elseif any(material==matOVC)
%                         OutPowerOVCfull(freqs)=OutPowerOVCfull(freqs)+0.5*condOVC*sigma(material)*abs((EF'*EF))*intw(pp)*det;
%                     end
%                 end

                if any(material==mat4K)
                    OutPower4Kfull=OutPower4Kfull+((0.5*cond4K*sigma(material)*intw(pp)*det).*EF_full);
                elseif any(material==mat77K)
                    OutPower77Kfull=OutPower77Kfull+((0.5*cond77K*sigma(material)*intw(pp)*det).*EF_full);
                elseif any(material==matOVC)
                    OutPowerOVCfull=OutPowerOVCfull+((0.5*condOVC*sigma(material)*intw(pp)*det).*EF_full);
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
                else
                gph(1:gesizet,1)=gphx(pp,1:gesizet)';
                gph(1:gesizet,2)=gphy(pp,1:gesizet)';
                gph(1:gesizet,3)=gphz(pp,1:gesizet)';
                
                % evaluate covairant mapping
                [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
                end
                
              % use stored functions
                ph(1:esizet,1:3)=(ephx(pp,1:esizet)'*axi(1:3))+...
                    (ephy(pp,1:esizet)'*aeta(1:3))+...
                    (ephz(pp,1:esizet)'*azeta(1:3));
                
                  cph(1:esizet,1:3)=(ecphx(pp,1:esizet)'*asxi(1:3))+...
                    (ecphy(pp,1:esizet)'*aseta(1:3))+...
                    (ecphz(pp,1:esizet)'*aszeta(1:3));
                
                
                % Compute the magnetic vector Potential
                efull = ph'*ldyn; % This is the magnetic vector potential
                
                
                % Compute the static magnetic flux density
                curlADC=cph'*lsolStatic;
            
                H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
                H1bas3D=zeros(3,3*esizeH1);

            
                H1bas3D(1,1:3:end)=H1bas;
                H1bas3D(2,2:3:end)=H1bas;
                H1bas3D(3,3:3:end)=H1bas;

            
                % compute the solution for this problem, at this integration point
                % in this element
                u_Dfull = H1bas3D*ldynmech;
                
                % Compute the electric field     

                curlACDrep = repmat(curlADC, 1, length(freqOut));
                crossfull = cross(curlACDrep,u_Dfull);
                ElectricFieldfull = (1i*(omegafull.')).*((crossfull)-(efull));
                EF_full = abs(diag(ElectricFieldfull'*ElectricFieldfull));

                % The output power is:
       
%                 for freqs = 1:length(freqOut)
%                   EF = ElectricFieldfull(:,freqs);
%                     if any(material==mat4K)
%                         OutPower4Kfull(freqs)=OutPower4Kfull(freqs)+0.5*cond4K*sigma(material)*abs((EF'*EF))*intw(pp)*det;
%                     elseif any(material==mat77K)
%                         OutPower77Kfull(freqs)=OutPower77Kfull(freqs)+0.5*cond77K*sigma(material)*abs((EF'*EF))*intw(pp)*det;
%                     elseif any(material==matOVC)
%                         OutPowerOVCfull(freqs)=OutPowerOVCfull(freqs)+0.5*condOVC*sigma(material)*abs((EF'*EF))*intw(pp)*det;
%                     end
%                 end

                                 
                if any(material==mat4K)
                    OutPower4Kfull=OutPower4Kfull+((0.5*cond4K*sigma(material)*intw(pp)*det).*EF_full);
                elseif any(material==mat77K)
                    OutPower77Kfull=OutPower77Kfull+((0.5*cond77K*sigma(material)*intw(pp)*det).*EF_full);
                elseif any(material==matOVC)
                    OutPowerOVCfull=OutPowerOVCfull+((0.5*condOVC*sigma(material)*intw(pp)*det).*EF_full);
                end

                
            end
            
        end
    end
    % end of loop over elements
end