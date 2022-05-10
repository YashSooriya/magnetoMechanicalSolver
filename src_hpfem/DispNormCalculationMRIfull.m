function [DispNorm4Kfull,DispNorm77Kfull,DispNormOVCfull] = DispNormCalculationMRIfull(mesh,unknown,Basis,Quadrature,probdata,Dynamic,freqOut)


% Extract relevant data from mesh structure

mat=mesh.mat;
esizeH1=probdata.esizeH1;
gorder=probdata.jb.gorder;

DispNorm4Kfull=zeros(length(freqOut),1);
DispNorm77Kfull=zeros(length(freqOut),1);
DispNormOVCfull=zeros(length(freqOut),1);

matOVC=probdata.matOVC;
mat77K=probdata.mat77K;
mat4K=probdata.mat4K;

%==================================================================================================
% Calculate Displacement norm
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
gphQuadx=Basis.gphQuadx;
gphQuady=Basis.gphQuady;
gphQuadz=Basis.gphQuadz;

% Extract relevant data from Quadrature structure
intw=Quadrature.intw;
nip=Quadrature.nip;


% Initialise variables
lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);



%--------------------------------------------------------------------------
% Loop over elements
%--------------------------------------------------------------------------
for i=1:nelem
    % Local to global mapping
    aaa=mapL2G_e(i);
    material=mat(aaa);
    
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
    ldyn = zeros(3*esizeH1, length(freqOut));
    
    for freqs = 1:length(freqOut)
            ldyn(:,freqs) = Dynamic(lunkv,freqs);
    end
    
    % choose correct set of stored basis functions
    if eltype(i)==1
        gphx=gphx1;
        gphy=gphy1;
        gphz=gphz1;
        phh1=phh11;
        H1basis=H1basis1;
    else
        gphx=gphx2;
        gphy=gphy2;
        gphz=gphz2;
        phh1=phh12;
        H1basis=H1basis2;
    end
    
    
    %-------------------------------------------------------------------------
    % evaluate covairant mapping (for linear geometry mapping is constant)
    if flag==0
        gph(1:gesizet,1)=gphx(1,1:gesizet)';
        gph(1:gesizet,2)=gphy(1,1:gesizet)';
        gph(1:gesizet,3)=gphz(1,1:gesizet)';
        
        [~,~,~,~,~,~,det]=jacobian_pre(flag,gesizet,gph,mycoord);
        
        
        for pp=1:nip
            
            ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
            
            % use stored functions
            
            H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
            H1bas3D=zeros(3,3*esizeH1);

            
            H1bas3D(1,1:3:end)=H1bas;
            H1bas3D(2,2:3:end)=H1bas;
            H1bas3D(3,3:3:end)=H1bas;

            
            % compute the solution for this problem, at this integration point
            % in this element
            u_Dfull = H1bas3D*ldyn;
               
                for freqs = 1:length(freqOut)
                    u_D = u_Dfull(:,freqs);
                    if any(material==mat4K)
                        DispNorm4Kfull(freqs)=DispNorm4Kfull(freqs)+intw(pp)*det*abs(u_D'*u_D);
                    elseif any(material==mat77K)
                        DispNorm77Kfull(freqs)=DispNorm77Kfull(freqs)+intw(pp)*det*abs(u_D'*u_D);
                    elseif any(material==matOVC)
                        DispNormOVCfull(freqs)=DispNormOVCfull(freqs)+intw(pp)*det*abs(u_D'*u_D);
                    end
                end
        
            
            
        end
        
    else
        
        for pp=1:nip
                if gorder==1
                gphQuad(1:gesizet,1)=gphQuadx(pp,1:gesizet)';
                gphQuad(1:gesizet,2)=gphQuady(pp,1:gesizet)';
                gphQuad(1:gesizet,3)=gphQuadz(pp,1:gesizet)';
                
                
                % evaluate covairant mapping
                [~,~,~,~,~,~,det]=jacobian_pre(flag,gesizet,gphQuad,mycoord);
                else
                gph(1:gesizet,1)=gphx(pp,1:gesizet)';
                gph(1:gesizet,2)=gphy(pp,1:gesizet)';
                gph(1:gesizet,3)=gphz(pp,1:gesizet)';
                ph1(1:gesizet,1)=phh1(pp,1:gesizet)';
                
                % evaluate covairant mapping
                [~,~,~,~,~,~,det]=jacobian_pre(flag,gesizet,gph,mycoord);
                end
            
            % use stored functions
            
            H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
            H1bas3D=zeros(3,3*esizeH1);

            
            H1bas3D(1,1:3:end)=H1bas;
            H1bas3D(2,2:3:end)=H1bas;
            H1bas3D(3,3:3:end)=H1bas;
            
          
            % compute the solution for this problem, at this integration point
            % in this element
            u_Dfull = H1bas3D*ldyn;
               
                for freqs = 1:length(freqOut)
                    u_D = u_Dfull(:,freqs);
                    if any(material==mat4K)
                        DispNorm4Kfull(freqs)=DispNorm4Kfull(freqs)+intw(pp)*det*abs(u_D'*u_D);
                    elseif any(material==mat77K)
                        DispNorm77Kfull(freqs)=DispNorm77Kfull(freqs)+intw(pp)*det*abs(u_D'*u_D);
                    elseif any(material==matOVC)
                        DispNormOVCfull(freqs)=DispNormOVCfull(freqs)+intw(pp)*det*abs(u_D'*u_D);
                    end
                end
                
            
            
        end
        
    end
    
    % end of loop over elements
end

DispNorm4Kfull=sqrt(DispNorm4Kfull);
DispNorm77Kfull=sqrt(DispNorm77Kfull);
DispNormOVCfull=sqrt(DispNormOVCfull);