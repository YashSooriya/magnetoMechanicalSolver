function [Matrices] =postProcessMatrices(Mesh,Unknown,Basis,Quadrature,ProblemData,Static,UnknownStatic)

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
nSubdoms=max(mat);
subFlag=Mesh.subFlag;

% Extract relevant data from ProblemData structure
esizet=ProblemData.esizet;
gorder=ProblemData.jb.gorder;
esizeH1=ProblemData.esizeH1;
sizeMech=3*esizeH1;
% Extract data from Static structure
solStatic=Static.sol(:,1);
nSparse=Static.nSparse;

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

% Extract relevant data from Unknown structure
lenEM=Unknown.EM.nunkt;
lenM=Unknown.system.nunkt-lenEM;


lec = zeros(6,gorder,3);
lfc = zeros(4,(gorder*(gorder-1)/2),3);



sigma=ProblemData.matr.sigma;
rho=ProblemData.matr.rho;

I_aa   = zeros(nSparse,nSubdoms);
J_aa   = zeros(nSparse,nSubdoms);
V_aa   = zeros(nSparse,nSubdoms);
I_uu   = zeros(nSparse,nSubdoms);
J_uu   = zeros(nSparse,nSubdoms);
V_uu   = zeros(nSparse,nSubdoms);
I_au   = zeros(nSparse,nSubdoms);
J_au   = zeros(nSparse,nSubdoms);
V_au   = zeros(nSparse,nSubdoms);
I_M    = zeros(nSparse,nSubdoms);
J_M    = zeros(nSparse,nSubdoms);
V_M    = zeros(nSparse,nSubdoms);
nz_aa  = zeros(1,nSubdoms);
nz_uu  = zeros(1,nSubdoms);
nz_au  = zeros(1,nSubdoms);
nz_M   = zeros(1,nSubdoms);




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
        lunkv2=Unknown.Mech.unknownsD(:,i)-lenEM-Unknown.EM.npec;
        lunkvStatic=UnknownStatic.EM.unknowns(:,i);
        
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
        
    P_AU=zeros(esizet,sizeMech);
    P_AA=zeros(esizet,esizet);
    P_UU=zeros(sizeMech,sizeMech);
    kE=zeros(sizeMech,sizeMech);
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
                
                
            
                H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
                H1bas3D=zeros(3,3*esizeH1);

            
                H1bas3D(1,1:3:end)=H1bas;
                H1bas3D(2,2:3:end)=H1bas;
                H1bas3D(3,3:3:end)=H1bas;

                % Compute BDC
                curlADC=cph'*lsolStatic;
                
        %==========================================================================
        % Vectorised implementation
        %==========================================================================
        B1=[0 curlADC(3) -curlADC(2); -curlADC(3) 0 curlADC(1); curlADC(2) -curlADC(1) 0];
        P_AU=P_AU+0.5*sigma(material)*ph*B1*H1bas3D*det*intw;
        P_AA=P_AA+0.5*sigma(material)*(ph*ph')*det*intw;
        P_UU=P_UU+0.5*sigma(material)*((curlADC*curlADC.')*(H1bas3D*H1bas3D')-(H1bas3D*curlADC.')*(H1bas3D*curlADC.').')*det*intw;
        %==========================================================================

                
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
                
                
            
                H1bas(1:esizeH1)=H1basis(pp,1:esizeH1);
            
                H1bas3D=zeros(3,3*esizeH1);

            
                H1bas3D(1,1:3:end)=H1bas;
                H1bas3D(2,2:3:end)=H1bas;
                H1bas3D(3,3:3:end)=H1bas;

                % Compute BDC
                curlADC=cph'*lsolStatic;
                
        %==========================================================================
        % Vectorised implementation
        %==========================================================================
        B1=[0 curlADC(3) -curlADC(2); -curlADC(3) 0 curlADC(1); curlADC(2) -curlADC(1) 0];
        P_AU=P_AU+0.5*sigma(material)*(ph*B1)*H1bas3D*det*intw(pp);
        P_AA=P_AA+0.5*sigma(material)*(ph*ph')*det*intw(pp);
        P_UU=P_UU+0.5*sigma(material)*((curlADC.'*curlADC)*(H1bas3D'*H1bas3D)-(H1bas3D'*curlADC)*(H1bas3D'*curlADC).')*det*intw(pp);
        kE=kE+0.5*rho(material)*(H1bas3D'*H1bas3D)*det*intw(pp);
        %==========================================================================


                
            end

        [I_aa(:,material),J_aa(:,material),V_aa(:,material),nz_aa(:,material)] = ...
            assembleMatrixComps(P_AA,I_aa(:,material),J_aa(:,material),V_aa(:,material),nz_aa(:,material),lunkv,lunkv,esizet,esizet);
          
        [I_uu(:,material),J_uu(:,material),V_uu(:,material),nz_uu(:,material)] = ...
            assembleMatrixComps(P_UU,I_uu(:,material),J_uu(:,material),V_uu(:,material),nz_uu(:,material),lunkv2,lunkv2,sizeMech,sizeMech);
        
        [I_au(:,material),J_au(:,material),V_au(:,material),nz_au(:,material)] = ...
            assembleMatrixComps(P_AU,I_au(:,material),J_au(:,material),V_au(:,material),nz_au(:,material),lunkv,lunkv2,esizet,sizeMech);
        
        [I_M(:,material),J_M(:,material),V_M(:,material),nz_M(:,material)]     = ...
            assembleMatrixComps(kE,I_M(:,material),J_M(:,material),V_M(:,material),nz_M(:,material),lunkv2,lunkv2,sizeMech,sizeMech);
            
        end
    end
    % end of loop over elements
end
Paa     = cell(1,nSubdoms);
Puu     = cell(1,nSubdoms);
Pau     = cell(1,nSubdoms);
kEnergy = cell(1,nSubdoms);


for i = 1:nSubdoms
    % Store the Power output matrix
    Paa{i} = sparse(I_aa(1:nz_aa(i),i),J_aa(1:nz_aa(i),i),V_aa(1:nz_aa(i),i),lenEM+Unknown.EM.npec,lenEM+Unknown.EM.npec);
    Puu{i} = sparse(I_uu(1:nz_uu(i),i),J_uu(1:nz_uu(i),i),V_uu(1:nz_uu(i),i),lenM+Unknown.Mech.npec,lenM+Unknown.Mech.npec);
    Pau{i} = sparse(I_au(1:nz_au(i),i),J_au(1:nz_au(i),i),V_au(1:nz_au(i),i),lenEM+Unknown.EM.npec,lenM+Unknown.Mech.npec);

    
    % Styore the kinetic energy matrix
    kEnergy{i} = sparse(I_M(1:nz_M(i),i),J_M(1:nz_M(i),i),V_M(1:nz_M(i),i),lenM+Unknown.Mech.npec,lenM+Unknown.Mech.npec);
end

Matrices.kEnergy  = kEnergy;
Matrices.Paa      = Paa;
Matrices.Puu      = Puu;
Matrices.Pau      = Pau;
Matrices.nSubdoms = nSubdoms;