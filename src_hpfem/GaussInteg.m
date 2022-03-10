function [M,C,K,Res]=GaussInteg(Basis,Quadrature,ProblemData,flag,mue,Adc,A0,mycoord,eltype,probstatic,material,subFlag,ielem,lec,lfc,xy,cond,coupling,FlagVolume,SourceMapping)
%=========================================================================
% Determine the elemental matrices
%=========================================================================
% Extract the relevant data from the input structures
if eltype==1
    ephx=Basis.ephx1;
    ephy=Basis.ephy1;
    ephz=Basis.ephz1;
    ecphx=Basis.ecphx1;
    ecphy=Basis.ecphy1;
    ecphz=Basis.ecphz1;
    gphx=Basis.gphx1;
    gphy=Basis.gphy1;
    gphz=Basis.gphz1;
    phh1=Basis.phh11;
    H1basis=Basis.H1basis1;
    gradH1basisx=Basis.gradH1basis1x;
    gradH1basisy=Basis.gradH1basis1y;
    gradH1basisz=Basis.gradH1basis1z;
    ephdfx=Basis.Face.ephdfx1;
    ephdfy=Basis.Face.ephdfy1;
    ephdfz=Basis.Face.ephdfz1;
    phh1f=Basis.Face.phh1f1;
    phH1f=Basis.Face.phH1f1;
    gphfx=Basis.Face.gphfx1;
    gphfy=Basis.Face.gphfy1;
    gphfz=Basis.Face.gphfz1;
elseif eltype==2
    ephx=Basis.ephx2;
    ephy=Basis.ephy2;
    ephz=Basis.ephz2;
    ecphx=Basis.ecphx2;
    ecphy=Basis.ecphy2;
    ecphz=Basis.ecphz2;
    gphx=Basis.gphx2;
    gphy=Basis.gphy2;
    gphz=Basis.gphz2;
    phh1=Basis.phh12;
    H1basis=Basis.H1basis2;
    gradH1basisx=Basis.gradH1basis2x;
    gradH1basisy=Basis.gradH1basis2y;
    gradH1basisz=Basis.gradH1basis2z;
    ephdfx=Basis.Face.ephdfx2;
    ephdfy=Basis.Face.ephdfy2;
    ephdfz=Basis.Face.ephdfz2;
    phh1f=Basis.Face.phh1f2;
    phH1f=Basis.Face.phH1f2;
    gphfx=Basis.Face.gphfx2;
    gphfy=Basis.Face.gphfy2;
    gphfz=Basis.Face.gphfz2;
end
gphQuadx=Basis.gphQuadx;
gphQuady=Basis.gphQuady;
gphQuadz=Basis.gphQuadz;
gphfQuadx=Basis.Face.gphQuadfx;
gphfQuady=Basis.Face.gphQuadfy;
gphfQuadz=Basis.Face.gphQuadfz;
pV=Basis.Face.p;
nmf=Basis.Face.nmf;
intxi=Quadrature.intxi;
inteta=Quadrature.inteta;
intzeta=Quadrature.intzeta;
intfw=Quadrature.intfw;
intfxi=Quadrature.intfxi;
intfet=Quadrature.intfet;
nipf=Quadrature.nipf;
esizet=ProblemData.esizet;
esizeH1=ProblemData.esizeH1;
gorder=ProblemData.jb.gorder;

%=========================================================================
% Initialisation of the linear system vectors/matrices
%=========================================================================

% Predefine the size of the blocks
R_A        = zeros(esizet,1);
R_U        = zeros(3*esizeH1,1);
% Predefine the size of the blocks
K_AA       = zeros(esizet);
K_AU       = zeros(esizet,3*esizeH1);
K_UA       = zeros(3*esizeH1,esizet);
K_UU       = zeros(3*esizeH1,3*esizeH1);
C_AA       = zeros(esizet);
C_AU       = zeros(esizet,3*esizeH1);
C_UA       = zeros(3*esizeH1,esizet);
C_UU       = zeros(3*esizeH1);
M_AA       = zeros(esizet);
M_AU       = zeros(esizet,3*esizeH1);
M_UA       = zeros(3*esizeH1,esizet);
M_UU       = zeros(3*esizeH1);



%--------------------------------------------------------------------------
% Define the elasticity tensor C_ijkl=D
%--------------------------------------------------------------------------

if subFlag==2
    D=ProblemData.matr.D{material};
else
    D=0;
end

%=========================================================================
% Compute elemental matrices volume contributions
%==========================================================================
gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;

if flag==0
    gesizet=4;
end

nip=Quadrature.nip;

% Initialize variables
B=zeros(6,3*esizeH1);
curlBasisEM=zeros(esizet,3);
BasisEM=zeros(esizet,3);
H1bas3D=zeros(3,3*esizeH1);
gradH1bas=zeros(esizeH1,3);
constantB=zeros(3*esizet,3);
A2=zeros(3*esizet,9);
Grad3D=zeros(3*esizeH1,9);
sigmae=zeros(9,1);
S2=zeros(9,esizet);
H1bas=zeros(1,esizeH1);
gph=zeros(gesizet,3);
bodyForce=zeros(3,1);
%========================================================================
% Compute the body terms
%========================================================================

if flag==0   % Linear geometry
    
    % Reorder necessary basis for jacobian calculation
    gph(1:gesizet,1)=gphx(1,1:gesizet)';
    gph(1:gesizet,2)=gphy(1,1:gesizet)';
    gph(1:gesizet,3)=gphz(1,1:gesizet)';
    
    % Compute jacobian
    [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
    
    for i=1:nip
        ph1(1:gesizet,1)=phh1(i,1:gesizet)';
        
        % Compute x,y,z
        [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
        
        % Compute the current source value
        srcfun=ProblemData.es.srcfun;
        arg=ProblemData.es.srcfunarg;
        current=srcfun(x,y,z,arg,material,probstatic);
        intw=Quadrature.intw(i);
        %curlADC2=[3*cos(3*x)*sin(4*y)*exp(5*z);4*sin(3*x)*cos(4*y)*exp(5*z);5*sin(3*x)*sin(4*y)*exp(5*z)];
        %----------------------------------------------------------------
        % Define Shape functions
        %----------------------------------------------------------------
        % Deine H(curl) basis
        BasisEM(1:esizet,1:3)=(ephx(i,1:esizet)'*axi(1:3))+...
            (ephy(i,1:esizet)'*aeta(1:3))+...
            (ephz(i,1:esizet)'*azeta(1:3));
        
        % Define curl of H(curl) basis
        curlBasisEM(1:esizet,1:3)=(ecphx(i,1:esizet)'*asxi(1:3))+...
            (ecphy(i,1:esizet)'*aseta(1:3))+...
            (ecphz(i,1:esizet)'*aszeta(1:3));
        
            %========================================================================      
            % Mechanical subdomain
            %========================================================================
             if FlagVolume==1
            
            % Compute source term for mechanical problem (body force)
            srcfun=ProblemData.es.srcfunmech;
            arg=ProblemData.es.srcfunargmech;
            bodyForce=srcfun(x,y,z,arg,material);
            
            % H1 basis functions
            H1bas(1:esizeH1)=H1basis(i,1:esizeH1);
            
            % Define H1 basis matrix for vectorised calculus
            H1bas3D(1,1:3:end)=H1bas;
            H1bas3D(2,2:3:end)=H1bas;
            H1bas3D(3,3:3:end)=H1bas;
            
            % Define gradient of H1 basis
            gradH1bas(1:esizeH1,1:3)=(gradH1basisx(i,1:esizeH1)'*axi(1:3))+...
                (gradH1basisy(i,1:esizeH1)'*aeta(1:3))+...
                (gradH1basisz(i,1:esizeH1)'*azeta(1:3));
            
            % Define matrix for vectorised calculus
            Grad3D(1:3:end,1)=gradH1bas(:,1);
            Grad3D(1:3:end,2)=gradH1bas(:,2);
            Grad3D(1:3:end,3)=gradH1bas(:,3);
            Grad3D(2:3:end,4)=gradH1bas(:,1);
            Grad3D(2:3:end,5)=gradH1bas(:,2);
            Grad3D(2:3:end,6)=gradH1bas(:,3);
            Grad3D(3:3:end,7)=gradH1bas(:,1);
            Grad3D(3:3:end,8)=gradH1bas(:,2);
            Grad3D(3:3:end,9)=gradH1bas(:,3);
            
            % Define B term for stiffnes calulation (Voigt notation)
            B(1,1:3:end)=gradH1bas(:,1);
            B(2,2:3:end)=gradH1bas(:,2);
            B(3,3:3:end)=gradH1bas(:,3);
            B(4,1:3:end)=gradH1bas(:,2);
            B(4,2:3:end)=gradH1bas(:,1);
            B(5,1:3:end)=gradH1bas(:,3);
            B(5,3:3:end)=gradH1bas(:,1);
            B(6,2:3:end)=gradH1bas(:,3);
            B(6,3:3:end)=gradH1bas(:,2);
        
             end
        %-----------------------------------------------------------------
        % Construct the directional derivative blocks
        %-----------------------------------------------------------------
      % Build each matrix (stiffness, damping, mass and residual)
        [K_AA,K_AU,K_UA,K_UU]=stiffLinear(K_AA,K_AU,K_UA,K_UU,esizet,mue,D,subFlag,Adc,det,curlBasisEM,B,intw,constantB,A2,Grad3D,S2,coupling,FlagVolume);
        [C_AA,C_AU]=dampLinear(C_AA,C_AU,esizet,det,BasisEM,intw,curlBasisEM,Adc,probstatic,subFlag,coupling,ProblemData,material,H1bas3D);
        [M_UU]=massLinear(M_UU,subFlag,ProblemData,material,det,H1bas3D,intw);
        [R_A,R_U]=residual(R_A,R_U,mue,Adc,A0,probstatic,det,BasisEM,curlBasisEM,H1bas3D,intw,current,bodyForce,Grad3D,sigmae,coupling,FlagVolume,SourceMapping);
        
    end
    
    %========================================================================================================================================
    % Compute the boundary terms (Residual and stiffness)
    %========================================================================================================================================
    
    % set up normals on the faces of the refeerence element
    ph1=zeros(gesizet,1);
    

%=========================================================================================
% Implement Neumann BC
%=========================================================================================
   for j=1:4
        if cond(ielem,j)==3
                     
            p=pV(:,:,j);
            
            % Initialise vectors
            Basis3Df=zeros(3*esizeH1,3);
            
            
            % Loop over integration points in face
            for pp=1:nipf
                
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
                [axi,aeta,azeta,~,~,~,~]=jacobian_pre(flag,gesizet,gph,mycoord) ;
                
                
                % Evaluate basis on faces
                
                Nhf=phH1f((j-1)*nipf+pp,1:esizeH1);
                
                ph(1:esizet,1:3)=(ephdfx((j-1)*nipf+pp,1:esizet)'*axi(1:3))+...
                    (ephdfy((j-1)*nipf+pp,1:esizet)'*aeta(1:3))+...
                    (ephdfz((j-1)*nipf+pp,1:esizet)'*azeta(1:3));
                
                
                
                % Create matrix for vectorised calculus
                Basis3Df(1:3:end,1)=Nhf;
                Basis3Df(2:3:end,2)=Nhf;
                Basis3Df(3:3:end,3)=Nhf;
                % computx x,y,z
                [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
                
                % Define normal vector
                nm(1)=(nmf(j,1)*axi(1))+(nmf(j,2)*aeta(1))+(nmf(j,3)*azeta(1));
                nm(2)=(nmf(j,1)*axi(2))+(nmf(j,2)*aeta(2))+(nmf(j,3)*azeta(2));
                nm(3)=(nmf(j,1)*axi(3))+(nmf(j,2)*aeta(3))+(nmf(j,3)*azeta(3));
                
                nd=sqrt((nm(1)^2)+(nm(2)^2)+(nm(3)^2));
                nm(1)=nm(1)/nd;
                nm(2)=nm(2)/nd;
                nm(3)=nm(3)/nd;
                
                
                % sign of computed normal is reversed from outward normal and so
                nm(1:3) = -nm(1:3);
                
                % Compute area for integration
                area=surdet(eltype,xi,eta,zeta,j,xy,lec,lfc,gorder+1);
                area=abs(area);
                % Compute curl(A) from Neumann function
                fun=ProblemData.es.neufun;
                arg=ProblemData.es.neufunarg;
                index=cond(ielem,j);
                curle=fun(x,y,z,index,arg,probstatic,ProblemData.matr.omega);
                
                if subFlag==2
                    % User problem file to define form of curl theta_i on the
                    % boundary/interface
                    fun=ProblemData.es.GradientExact;
                    arg=ProblemData.es.gradfunarg;
                    grad_exact=fun(x,y,z,arg);
                else
                    grad_exact=zeros(3,3);
                end
                
                
                % Compute the Neumann BC
                [R_A,R_U]= NeumannBC(R_A,R_U,intfw,mue,curle,esizet,nm,area,subFlag,ProblemData,grad_exact,ph,Basis3Df,material,pp);
            end
        end
    end
    
    
    
else   % High order geometry
    for i=1:nip
        if gorder>1
        gph(1:gesizet,1)=gphx(i,1:gesizet)';
        gph(1:gesizet,2)=gphy(i,1:gesizet)';
        gph(1:gesizet,3)=gphz(i,1:gesizet)';
        
        % Mapping
        [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
        ph1(1:gesizet,1)=phh1(i,1:gesizet)';
        
        % Get physical coordinates
        [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
        else
        gphQuad(1:gesizet,1)=gphQuadx(i,1:10)';
        gphQuad(1:gesizet,2)=gphQuady(i,1:10)';
        gphQuad(1:gesizet,3)=gphQuadz(i,1:10)';
         % Mapping
        [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gphQuad,mycoord);
        % Get physical coordinates
        [x,y,z]= getxyzq(mycoord,intxi(i),inteta(i),intzeta(i));
        end
            
        %curlADC2=[3*cos(3*x)*sin(4*y)*exp(5*z);4*sin(3*x)*cos(4*y)*exp(5*z);5*sin(3*x)*sin(4*y)*exp(5*z)];
        % Use problem file to define src function for EM and Mech
        srcfun=ProblemData.es.srcfun;
        arg=ProblemData.es.srcfunarg;
        current=srcfun(x,y,z,arg,material,probstatic);
        intw=Quadrature.intw(i);
        
        % Define H(curl) basis
        BasisEM(1:esizet,1:3)=(ephx(i,1:esizet)'*axi(1:3))+...
            (ephy(i,1:esizet)'*aeta(1:3))+...
            (ephz(i,1:esizet)'*azeta(1:3));
        
        % Define curl of H(curl) basis
        curlBasisEM(1:esizet,1:3)=(ecphx(i,1:esizet)'*asxi(1:3))+...
            (ecphy(i,1:esizet)'*aseta(1:3))+...
            (ecphz(i,1:esizet)'*aszeta(1:3));
        
        
        
        
            % Mechanical subdomain
            if FlagVolume==1
            % Define source term for mechanical problem (body force)
            srcfun=ProblemData.es.srcfunmech;
            arg=ProblemData.es.srcfunargmech;
            bodyForce=srcfun(x,y,z,arg,material);
            
            % Define H1 basis
            H1bas(1:esizeH1)=H1basis(i,1:esizeH1);
            
            % Define matrix for vectorised calculus
            H1bas3D(1,1:3:end)=H1bas;
            H1bas3D(2,2:3:end)=H1bas;
            H1bas3D(3,3:3:end)=H1bas;
            
            % Define gradient of H1 basis
            gradH1bas(1:esizeH1,1:3)=(gradH1basisx(i,1:esizeH1)'*axi(1:3))+...
                (gradH1basisy(i,1:esizeH1)'*aeta(1:3))+...
                (gradH1basisz(i,1:esizeH1)'*azeta(1:3));
            
            % Define matrix for vectorised calculus
            Grad3D(1:3:end,1)=gradH1bas(:,1);
            Grad3D(1:3:end,2)=gradH1bas(:,2);
            Grad3D(1:3:end,3)=gradH1bas(:,3);
            Grad3D(2:3:end,4)=gradH1bas(:,1);
            Grad3D(2:3:end,5)=gradH1bas(:,2);
            Grad3D(2:3:end,6)=gradH1bas(:,3);
            Grad3D(3:3:end,7)=gradH1bas(:,1);
            Grad3D(3:3:end,8)=gradH1bas(:,2);
            Grad3D(3:3:end,9)=gradH1bas(:,3);
            
            % Define B matrix (Voigt notation)
            B(1,1:3:end)=gradH1bas(:,1);
            B(2,2:3:end)=gradH1bas(:,2);
            B(3,3:3:end)=gradH1bas(:,3);
            B(4,1:3:end)=gradH1bas(:,2);
            B(4,2:3:end)=gradH1bas(:,1);
            B(5,1:3:end)=gradH1bas(:,3);
            B(5,3:3:end)=gradH1bas(:,1);
            B(6,2:3:end)=gradH1bas(:,3);
            B(6,3:3:end)=gradH1bas(:,2);
        
            end
        
        % Build each matrix (stiffness, damping, mass and residual)
        [K_AA,K_AU,K_UA,K_UU]=stiffLinear(K_AA,K_AU,K_UA,K_UU,esizet,mue,D,subFlag,Adc,det,curlBasisEM,B,intw,constantB,A2,Grad3D,S2,coupling,FlagVolume);
        [C_AA,C_AU]=dampLinear(C_AA,C_AU,esizet,det,BasisEM,intw,curlBasisEM,Adc,probstatic,subFlag,coupling,ProblemData,material,H1bas3D);
        [M_UU]=massLinear(M_UU,subFlag,ProblemData,material,det,H1bas3D,intw);
        [R_A,R_U]=residual(R_A,R_U,mue,Adc,A0,probstatic,det,BasisEM,curlBasisEM,H1bas3D,intw,current,bodyForce,Grad3D,sigmae,coupling,FlagVolume,SourceMapping);
        
    end

     
    % set up normals on the faces of the refeerence element
    
    ph1=zeros(gesizet,1);
         
%=========================================================================================
% Implement Neumann BC
%=========================================================================================
   for j=1:4
        if cond(ielem,j)==3
            
            
            
            p=pV(:,:,j);
            
            % Initialise vectors
            Basis3Df=zeros(3*esizeH1,3);
            
            
            % Loop over integration points in face
            for pp=1:nipf
                
                                l(1)=1-intfxi(pp)-intfet(pp);
                l(2)=intfxi(pp);
                l(3)=intfet(pp);
                
                xi = l(1:3)*p(1:3,1);
                eta = l(1:3)*p(1:3,2);
                zeta = l(1:3)*p(1:3,3);
                
         if gorder>1
        gph(1:gesizet,1)=gphfx(pp,1:gesizet)';
        gph(1:gesizet,2)=gphfy(pp,1:gesizet)';
        gph(1:gesizet,3)=gphfz(pp,1:gesizet)';
        
        % Mapping
        [axi,aeta,azeta,~,~,~,~]=jacobian_pre(flag,gesizet,gph,mycoord);
        ph1(1:gesizet,1)=phh1(i,1:gesizet)';
        
        % Get physical coordinates
        [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
        else
        gphQuad(1:gesizet,1)=gphfQuadx(pp,1:10)';
        gphQuad(1:gesizet,2)=gphfQuady(pp,1:10)';
        gphQuad(1:gesizet,3)=gphfQuadz(pp,1:10)';
         % Mapping
        [axi,aeta,azeta,~,~,~,~]=jacobian_pre(flag,gesizet,gphQuad,mycoord);
        % Get physical coordinates
        [x,y,z]= getxyzq(mycoord,xi,eta,zeta);
        end
        

              
                
                
                % Evaluate basis on faces
                
                Nhf=phH1f((j-1)*nipf+pp,1:esizeH1);
                
                ph(1:esizet,1:3)=(ephdfx((j-1)*nipf+pp,1:esizet)'*axi(1:3))+...
                    (ephdfy((j-1)*nipf+pp,1:esizet)'*aeta(1:3))+...
                    (ephdfz((j-1)*nipf+pp,1:esizet)'*azeta(1:3));
                
                
                
                % Create matrix for vectorised calculus
                Basis3Df(1:3:end,1)=Nhf;
                Basis3Df(2:3:end,2)=Nhf;
                Basis3Df(3:3:end,3)=Nhf;
           
                
                % Define normal vector
                nm(1)=(nmf(j,1)*axi(1))+(nmf(j,2)*aeta(1))+(nmf(j,3)*azeta(1));
                nm(2)=(nmf(j,1)*axi(2))+(nmf(j,2)*aeta(2))+(nmf(j,3)*azeta(2));
                nm(3)=(nmf(j,1)*axi(3))+(nmf(j,2)*aeta(3))+(nmf(j,3)*azeta(3));
                
                nd=sqrt((nm(1)^2)+(nm(2)^2)+(nm(3)^2));
                nm(1)=nm(1)/nd;
                nm(2)=nm(2)/nd;
                nm(3)=nm(3)/nd;
                
                
                % sign of computed normal is reversed from outward normal and so
                nm(1:3) = -nm(1:3);
                
                % Compute area for integration
                if gorder==1
                    area=surdetq(eltype,xi,eta,zeta,j,mycoord);
                else
                    area=surdet(eltype,xi,eta,zeta,j,xy,lec,lfc,gorder+1);
                end
                area=abs(area);
                % Compute curl(A) from Neumann function
                fun=ProblemData.es.neufun;
                arg=ProblemData.es.neufunarg;
                index=cond(ielem,j);
                curle=fun(x,y,z,index,arg,probstatic,ProblemData.matr.omega);
                
                if subFlag==2
                    % User problem file to define form of curl theta_i on the
                    % boundary/interface
                    fun=ProblemData.es.GradientExact;
                    arg=ProblemData.es.gradfunarg;
                    grad_exact=fun(x,y,z,arg);
                else
                    grad_exact=zeros(3,3);
                end
                
                
                % Compute the Neumann BC
                [R_A,R_U]= NeumannBC(R_A,R_U,intfw,mue,curle,esizet,nm,area,subFlag,ProblemData,grad_exact,ph,Basis3Df,material,pp);
            end
        end
    end
    
end


%=========================================================================
% Construct the local (element-wise) Linear system
%=========================================================================

% Construct the linearised mass matrix
M    = [M_AA M_AU;M_UA M_UU];

% Construct the linearised damping matrix
C    = [C_AA C_AU;C_UA C_UU];

% Construct the linearised stiffness matrix
K    = [K_AA K_AU;K_UA K_UU];

% Local residual vector
Res  = [R_A;R_U];

end