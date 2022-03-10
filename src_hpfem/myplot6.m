function myplot6(mesh,unknown,probdata,Solution,problem,probstatic,freq)


% Extract relevant data from mesh structure
coord=mesh.Coordinates;
intma=mesh.intma;
glob=mesh.edge.glob;
globfa=mesh.face.globfa;
eltype=mesh.eltype;
mat=mesh.mat;
edgecof=mesh.edgecof;
facecof=mesh.facecof;



order=probdata.order;
esizet=probdata.esizet;
gorder=probdata.jb.gorder;

% plot the result on a line through the domain
N=probdata.mesh.N;
point=probdata.mesh.point4;

sol=Solution.sol;
% interpolation__________________________________________________________________________

 % number of interpolation points



nelem=size(intma,1);
omega=2*pi*freq;
mu0=probdata.matr.muz;


cores = zeros(N,2);  % coresponding matrix
for i = 1:N
    in = -1;
    j=0;
    while in <0 && j<nelem
        %for j = 1:nelemt
        j=j+1;
        xyt=coord(intma(j,:),:);
        
        [in,onf]= intetra(xyt(1,:),xyt(2,:),xyt(3,:),xyt(4,:),point(i,:));
        if in ==1
            cores(i,:) = [i, j];            % cores = [point  element];
        end
    end
end

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

nrhs=probdata.jb.nrhs;
mypoint=[];
outh=zeros(N,3,nrhs);
fnd=0;
for i = 1:N
    if cores(i,2)~=0
        fnd=fnd+1;
        xyz=coord(intma(cores(i,2),1:4),1:3);
        
        Axy = [(xyz(2,1)-xyz(1,1))/2, (2*xyz(3,1)-xyz(2,1)-xyz(1,1))/(2*sqrt(3)), (3*xyz(4,1)-xyz(3,1)-xyz(2,1)-xyz(1,1))/(2*sqrt(6));
            (xyz(2,2)-xyz(1,2))/2, (2*xyz(3,2)-xyz(2,2)-xyz(1,2))/(2*sqrt(3)),  (3*xyz(4,2)-xyz(3,2)-xyz(2,2)-xyz(1,2))/(2*sqrt(6));
            (xyz(2,3)-xyz(1,3))/2, (2*xyz(3,3)-xyz(2,3)-xyz(1,3))/(2*sqrt(3)),  (3*xyz(4,3)-xyz(3,3)-xyz(2,3)-xyz(1,3))/(2*sqrt(6))];
        fxy = [point(i,1)-(xyz(1,1)+xyz(2,1))/2;
            point(i,2)-(xyz(1,2)+xyz(2,2))/2;
            point(i,3)-(xyz(1,3)+xyz(2,3))/2];
        pxieta = Axy\fxy;
        
        % local coordinate corresponding to plotting point
        xi = pxieta(1);
        eta = pxieta(2);
        zeta =pxieta(3);
        
        % check that this is a valid point inside the reference tet
        [in,onf]= intetra(v(1,:),v(2,:),v(3,:),v(4,:),pxieta');
        if in==1
            %disp('inside tetrahedron')
            mypoint=[mypoint;point(i,:)];
            
            %   compute mappings for these elements too
            flag=0 ;
            lec=zeros(6,gorder*3);
            lfc=zeros(4,3*gorder*(gorder-1)/2);
            for j=1:6
                for p=1:gorder
                    for k=1:3
                        lec(j,p,k)=edgecof(glob(cores(i,2),j),((p-1)*3)+k);
                        if abs(edgecof(glob(cores(i,2),j),((p-1)*3)+k))>0.00000001
                            flag=1;
                        end
                    end
                end
            end
            
            for j=1:4
                for p=1:gorder*(gorder-1)/2
                    for k=1:3
                        lfc(j,p,k)=facecof(globfa(cores(i,2),j),((p-1)*3)+k);
                    end
                end
            end
            
            
            gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
            
            mycoord = zeros(gesizet,3);
            
            %transfer coefficents to locations vertices
            mycoord(1:4,1:3) = xyz(1:4,1:3);
            
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
            

            % obtain mapping
            [axi,aeta,azeta,asxi,aseta,aszeta,~]=jacobian(xyz,xi,eta,zeta,...
                eltype(cores(i,2)),lec,lfc,flag,gorder,mycoord) ;
            
            % as point different for every element, no benefit in storing!
            if eltype(cores(i,2))==1
                phmxyz= curlbasis(xi,eta,zeta,asxi,aseta,aszeta,order,1,esizet);
                phmxyz1=basis(xi,eta,zeta,axi,aeta,azeta,order,1,esizet);
            else
                phmxyz= curlbasis(xi,eta,zeta,asxi,aseta,aszeta,order,2,esizet);
                 phmxyz1=basis(xi,eta,zeta,axi,aeta,azeta,order,2,esizet);
            end
            

                lunkv=unknown.EM.unknowns(:,cores(i,2));
            
            %     Magnetic Field
            he=zeros(esizet,1);
            he(lunkv>0) = sol(lunkv(lunkv>0),1);
            
            htil= phmxyz'*he/probdata.matr.mu(mat(cores(i,2)))/probdata.matr.muz;
            outh(fnd,:,:)=htil;
            sigma=probdata.matr.sigma(mat(cores(i,2)));
            Eddytil=phmxyz1'*he*(complex(0,-1)*omega*sigma);
            
            % It would not be difficult to extend this to compute E and J
            % at desired plotting points as well....
            
        end
        
    end
end

% plot out the solution
N=size(mypoint,1);


    
    outhrel=zeros(N,nrhs);
    for i = 1:N
        xp = mypoint(i,1);
        yp = mypoint(i,2);
        zp = mypoint(i,3);
        zcoord(i)=zp;
        r(i)=norm([xp,yp,zp]);
        x(i)=xp;
        problem=0;
        %   compute total B = mu H = curl A
        if problem==1 || problem==8
        fun=probdata.es.exactcurlfun;
        arg=probdata.es.exactcurlfunarg;
        [~,hf,H0]=fun(xp,yp,zp,arg,probstatic,omega);
        end
        %    store normalised perturbed H
        for j=1:nrhs
            if problem==1 
            outhrel(i,j)=norm(outh(i,:,j)*mu0-(H0(:,j).')*mu0)/(norm(H0(:,j)*mu0));
            hptrel(i,j)= norm(hf(:,j)*mu0-H0(:,j)*mu0)/(norm(H0(:,j)*mu0)); 
            bplo(i,j)=norm(outh(i,:,j)*mu0);
           
            elseif problem==8
            outhrel(i,j)=norm(outh(i,:,j)*mu0-[0;0;1].');
            hptrel(i,j)= norm(hf(:,j)*mu0-[0;0;1]); 
            bplo(i,j)=norm(outh(i,:,j)*mu0);
           
            else
            outhrel(i,j)=norm(outh(i,:,j)*mu0-[1;0;0].');
            bplo(i,j)=norm(outh(i,:,j)*mu0);
            bploz(i,j)=norm(outh(i,3,j)*mu0);
            end
            
            
            
            
        end
    end
   
   % Saving variables to plot for different p for TEAM 7
    
 
    
      
      
    
       
      
       
        
        figure
        plot(zcoord,bploz,'b-');
        xlabel('distance');
        ylabel('Magnetic field');
       
        save('LineOffAxis_p2Fine','zcoord','bploz');
        

        
  
     
        
     
     
        
    
    
    % Here we can add a possibility of just plotting out H, J or E as
    % desired.

