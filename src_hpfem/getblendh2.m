% Function builds coefficents for geometry approximation using an L2
% projection from an exact sphere geometry to a set of heirarchic polynomials of
% degree gorder+1

function [edgecof,facecof,mesh]=getblendh2(mesh,gorder,probdata)

% edgecof = blending coefficents for edge correction
% facecof = blending coefficents for face correction
% ae = local matrix used in determining edge coefficents
% se = local matrix used in determining edge coefficents
% re = local matrix used in determining edge coefficents
% af = local matrix used in determining face coefficents
% sf = local matrix used in determining face coefficents
% rf = local matrix used in determining face coefficents
% phef = vector of H^1 basis functions for geometry correction


% Extract relevant data from mesh structure
nelem=mesh.Nelements;
nside=mesh.edge.nside;
nface=mesh.face.nface;
intma=mesh.intma;
coord=mesh.Coordinates;
cond=mesh.face.cond;
glob=mesh.edge.glob;
globfa=mesh.face.globfa;
eltype=mesh.eltype;
mat=mesh.mat;


lec = zeros(6,gorder,3);
lfc = zeros(4,gorder*(gorder-1)/2,3);
in = 0;

gesizet = (gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
order = gorder;

% we need more than just (p+1)^2 as the ture curve is nonpoly
nipe = floor(max((2*(order+1)+1)/2,2));
% nipe = order*3

if nipe > 100
    error(message('increase dimension for int in getblend.m'));
end
maxnip = 5000;

% compute integration points on the interval (-1,1)
% call gauleg(-1,1,xie,we,nipe)
% int_-1^+1 dx

kind = 1;
kpts =0;
endpts(1) = -1;
endpts(2) = 1;

[~,xie,we]=gaussq1(kind,nipe,0,0,kpts,endpts);

% compute the intgration points over the face of an element
[intfxi,intfeta,intfw,nipf]=intpointf(nipe,maxnip);
% compute the intgration points over the volume of an element
[intxi,inteta,intzeta,intw,nipv]=intpoints(nipe,maxnip);

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
%store geometry basis functions over volume
gphx1=zeros(nipv,gesizet);
gphy1=zeros(nipv,gesizet);
gphz1=zeros(nipv,gesizet);

gphx2=zeros(nipv,gesizet);
gphy2=zeros(nipv,gesizet);
gphz2=zeros(nipv,gesizet);
for i=1:nipv
    gph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,gorder+1,1,gesizet);
    
    gphx1(i,1:gesizet)=gph(1:gesizet,1)';
    gphy1(i,1:gesizet)=gph(1:gesizet,2)';
    gphz1(i,1:gesizet)=gph(1:gesizet,3)';
    
    gph=gbasish1(intxi(i),inteta(i),intzeta(i),axi,aeta,azeta,gorder+1,2,gesizet);
    gphx2(i,1:gesizet)=gph(1:gesizet,1)';
    gphy2(i,1:gesizet)=gph(1:gesizet,2)';
    gphz2(i,1:gesizet)=gph(1:gesizet,3)';
end

safety = 1;
% vertices of reference element
v(1,1) = -1;
v(1,2) = 0;
v(1,3) = 0;

v(2,1) = 1;
v(2,2) = 0;
v(2,3) = 0;

v(3,1) = 0;
v(3,2) = sqrt(3);
v(3,3) = 0;

v(4,1) = 0;
v(4,2) = sqrt(3)/3;
v(4,3) = 2*(sqrt(2)/sqrt(3));

% zero edge coefficents
edgecof = zeros(nside,3*(order-1)+3);

% zero face coefficents
facecof = zeros(nface,(order*(order-1)/2-1)*3+3);

if gorder ~=0
    while 1
        rin = probdata.jb.rin(in+1);% str2num(jobdata{17+in*3});
        bcno = probdata.jb.sufv(in+1);%str2num(jobdata{18+in*3});
        % loop over elements in mesh
        for i =1:nelem % elem
            % local face to edge connections
            % lfedge contains the edge numbers on each face
            if eltype(i)==1
                lfedge(1,1)=1;
                lfedge(1,2)=5;
                lfedge(1,3)=6;
            else
                lfedge(1,1)=1;
                lfedge(1,2)=6;
                lfedge(1,3)=5;
            end
            lfedge(2,1)=2;
            lfedge(2,2)=4;
            lfedge(2,3)=6;
            
            lfedge(3,1)=3;
            lfedge(3,2)=4;
            lfedge(3,3)=5;
            
            if eltype(i)==1
                lfedge(4,1)=3;
                lfedge(4,2)=2;
                lfedge(4,3)=1;
            else
                lfedge(4,1)=2;
                lfedge(4,2)=3;
                lfedge(4,3)=1;
            end
            
            % loop over face in mesh
            for j = 1:4 %fac
                if cond(i,j)==bcno && mat(i)==2
                    % this is a boundary face and a blending correction should be applied
                    if cond(i,j) ==bcno
                        rs = rin;
                    else
                        rs = rout;
                    end
                    
                    % set up vertices of the face and move back to sphere surface
                    if j==1
                        p(1,1:3) = v(2,1:3);
                        p(2,1:3) = v(3,1:3);
                        p(3,1:3) = v(4,1:3);
                        
                        rd = sqrt((coord(intma(i,2),1)^2)+(coord(intma(i,2),2)^2)+...
                            (coord(intma(i,2),3)^2));
                        coord(intma(i,2),1)=(coord(intma(i,2),1)/rd)*rs;
                        coord(intma(i,2),2)=(coord(intma(i,2),2)/rd)*rs;
                        coord(intma(i,2),3)=(coord(intma(i,2),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,3),1)^2)+(coord(intma(i,3),2)^2)+...
                            (coord(intma(i,3),3)^2));
                        coord(intma(i,3),1)=(coord(intma(i,3),1)/rd)*rs;
                        coord(intma(i,3),2)=(coord(intma(i,3),2)/rd)*rs;
                        coord(intma(i,3),3)=(coord(intma(i,3),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,4),1)^2)+(coord(intma(i,4),2)^2)+...
                            (coord(intma(i,4),3)^2));
                        coord(intma(i,4),1)=(coord(intma(i,4),1)/rd)*rs;
                        coord(intma(i,4),2)=(coord(intma(i,4),2)/rd)*rs;
                        coord(intma(i,4),3)=(coord(intma(i,4),3)/rd)*rs;
                    elseif j==2
                        p(1,1:3) = v(3,1:3);
                        p(2,1:3) = v(1,1:3);
                        p(3,1:3) = v(4,1:3);
                        
                        rd = sqrt((coord(intma(i,1),1)^2)+(coord(intma(i,1),2)^2)+...
                            (coord(intma(i,1),3)^2));
                        coord(intma(i,1),1)=(coord(intma(i,1),1)/rd)*rs;
                        coord(intma(i,1),2)=(coord(intma(i,1),2)/rd)*rs;
                        coord(intma(i,1),3)=(coord(intma(i,1),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,3),1)^2)+(coord(intma(i,3),2)^2)+...
                            (coord(intma(i,3),3)^2));
                        coord(intma(i,3),1)=(coord(intma(i,3),1)/rd)*rs;
                        coord(intma(i,3),2)=(coord(intma(i,3),2)/rd)*rs;
                        coord(intma(i,3),3)=(coord(intma(i,3),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,4),1)^2)+(coord(intma(i,4),2)^2)+...
                            (coord(intma(i,4),3)^2));
                        coord(intma(i,4),1)=(coord(intma(i,4),1)/rd)*rs;
                        coord(intma(i,4),2)=(coord(intma(i,4),2)/rd)*rs;
                        coord(intma(i,4),3)=(coord(intma(i,4),3)/rd)*rs;
                    elseif j==3
                        p(1,1:3) = v(1,1:3);
                        p(2,1:3) = v(2,1:3);
                        p(3,1:3) = v(4,1:3);
                        
                        rd = sqrt((coord(intma(i,1),1)^2)+(coord(intma(i,1),2)^2)+...
                            (coord(intma(i,1),3)^2));
                        coord(intma(i,1),1)=(coord(intma(i,1),1)/rd)*rs;
                        coord(intma(i,1),2)=(coord(intma(i,1),2)/rd)*rs;
                        coord(intma(i,1),3)=(coord(intma(i,1),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,2),1)^2)+(coord(intma(i,2),2)^2)+...
                            (coord(intma(i,2),3)^2));
                        coord(intma(i,2),1)=(coord(intma(i,2),1)/rd)*rs;
                        coord(intma(i,2),2)=(coord(intma(i,2),2)/rd)*rs;
                        coord(intma(i,2),3)=(coord(intma(i,2),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,4),1)^2)+(coord(intma(i,4),2)^2)+...
                            (coord(intma(i,4),3)^2));
                        coord(intma(i,4),1)=(coord(intma(i,4),1)/rd)*rs;
                        coord(intma(i,4),2)=(coord(intma(i,4),2)/rd)*rs;
                        coord(intma(i,4),3)=(coord(intma(i,4),3)/rd)*rs;
                    else
                        p(1,1:3) = v(1,1:3);
                        p(2,1:3) = v(3,1:3);
                        p(3,1:3) = v(2,1:3);
                        
                        rd = sqrt((coord(intma(i,2),1)^2)+(coord(intma(i,2),2)^2)+...
                            (coord(intma(i,2),3)^2));
                        coord(intma(i,2),1)=(coord(intma(i,2),1)/rd)*rs;
                        coord(intma(i,2),2)=(coord(intma(i,2),2)/rd)*rs;
                        coord(intma(i,2),3)=(coord(intma(i,2),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,3),1)^2)+(coord(intma(i,3),2)^2)+...
                            (coord(intma(i,3),3)^2));
                        coord(intma(i,3),1)=(coord(intma(i,3),1)/rd)*rs;
                        coord(intma(i,3),2)=(coord(intma(i,3),2)/rd)*rs;
                        coord(intma(i,3),3)=(coord(intma(i,3),3)/rd)*rs;
                        
                        rd = sqrt((coord(intma(i,1),1)^2)+(coord(intma(i,1),2)^2)+...
                            (coord(intma(i,1),3)^2));
                        coord(intma(i,1),1)=(coord(intma(i,1),1)/rd)*rs;
                        coord(intma(i,1),2)=(coord(intma(i,1),2)/rd)*rs;
                        coord(intma(i,1),3)=(coord(intma(i,1),3)/rd)*rs;
                    end
                    
                    % transfer coordinates to local array
                    xy = coord(intma(i,1:4),1:3);
                    
                    % loop over integration points on an edge
                    for iedge = 1:3 % ed
                        % zero matrices for accumalation
                        a = zeros(order,order);
                        r = zeros(order,3);
                        
                        % loop over integration points on an edge
                        for pp = 1:nipe % ip
                            if lfedge(j,iedge)==1
                                % xit,etat, zetat is a point on the edge
                                if eltype(i) ==1
                                    xit = 0.5*((v(2,1)*(1-xie(pp)))+(v(3,1)*(1+xie(pp))));
                                    etat = 0.5*((v(2,2)*(1-xie(pp)))+(v(3,2)*(1+xie(pp))));
                                    zetat = 0.5*((v(2,3)*(1-xie(pp)))+(v(3,3)*(1+xie(pp))));
                                else
                                    xit = 0.5*((v(3,1)*(1-xie(pp)))+(v(2,1)*(1+xie(pp))));
                                    etat = 0.5*((v(3,2)*(1-xie(pp)))+(v(2,2)*(1+xie(pp))));
                                    zetat = 0.5*((v(3,3)*(1-xie(pp)))+(v(2,3)*(1+xie(pp))));
                                end
                            elseif lfedge(j,iedge)==2
                                xit = 0.5*((v(1,1)*(1-xie(pp)))+(v(3,1)*(1+xie(pp))));
                                etat = 0.5*((v(1,2)*(1-xie(pp)))+(v(3,2)*(1+xie(pp))));
                                zetat = 0.5*((v(1,3)*(1-xie(pp)))+(v(3,3)*(1+xie(pp))));
                            elseif lfedge(j,iedge)==3
                                xit = 0.5*((v(1,1)*(1-xie(pp)))+(v(2,1)*(1+xie(pp))));
                                etat = 0.5*((v(1,2)*(1-xie(pp)))+(v(2,2)*(1+xie(pp))));
                                zetat = 0.5*((v(1,3)*(1-xie(pp)))+(v(2,3)*(1+xie(pp))));
                            elseif lfedge(j,iedge)==4
                                xit = 0.5*((v(1,1)*(1-xie(pp)))+(v(4,1)*(1+xie(pp))));
                                etat = 0.5*((v(1,2)*(1-xie(pp)))+(v(4,2)*(1+xie(pp))));
                                zetat = 0.5*((v(1,3)*(1-xie(pp)))+(v(4,3)*(1+xie(pp))));
                            elseif lfedge(j,iedge)==5
                                xit = 0.5*((v(2,1)*(1-xie(pp)))+(v(4,1)*(1+xie(pp))));
                                etat = 0.5*((v(2,2)*(1-xie(pp)))+(v(4,2)*(1+xie(pp))));
                                zetat = 0.5*((v(2,3)*(1-xie(pp)))+(v(4,3)*(1+xie(pp))));
                            elseif lfedge(j,iedge)==6
                                xit = 0.5*((v(3,1)*(1-xie(pp)))+(v(4,1)*(1+xie(pp))));
                                etat = 0.5*((v(3,2)*(1-xie(pp)))+(v(4,2)*(1+xie(pp))));
                                zetat = 0.5*((v(3,3)*(1-xie(pp)))+(v(4,3)*(1+xie(pp))));
                            end
                            
                            %compute H1 basis for this edge
                            ph=basish1(gesizet,xit,etat,zetat,gorder+1,eltype(i));
                            
                            % initialize correction term
                            gc = ph(1:4,1)'*xy(1:4,1:3);
                            
                            % compute the correct coordinates on the surface of the sphere
                            rd=sqrt((gc(1)^2)+(gc(2)^2)+(gc(3)^2));
                            x=gc(1)/rd*rin;
                            y=gc(2)/rd*rin;
                            z=gc(3)/rd*rin;
                            for ii=1:order
                                row=4+lfedge(j,iedge)+6*(ii-1);
                                for jj=1:order
                                    col=4+lfedge(j,iedge)+6*(jj-1);
                                    a(ii,jj)=a(ii,jj)+(we(pp)*ph(row)*ph(col));
                                end
                                r(ii,1)=r(ii,1)+(ph(row)*we(pp)*(x-gc(1)));
                                r(ii,2)=r(ii,2)+(ph(row)*we(pp)*(y-gc(2)));
                                r(ii,3)=r(ii,3)+(ph(row)*we(pp)*(z-gc(3)));
                            end
                        end %ip
                        % solve the system for the g coefficents
                        sol = a\r;
                        
                        for jj=1:order
                            for k=1:3
                                if edgecof(glob(i,lfedge(j,iedge)),((jj-1)*3)+k)~=0 && ...
                                        abs(edgecof(glob(i,lfedge(j,iedge)),((jj-1)*3)+k)-...
                                        sol(jj,k)*safety)>0.0001
                                    if abs(edgecof(glob(i,lfedge(j,iedge)),((jj-1)*3)+k))>...
                                            abs(sol(jj,k)*safety)
                                        edgecof(glob(i,lfedge(j,iedge)),((jj-1)*3)+k)= sol(jj,k)*safety;
                                    end
                                else
                                    edgecof(glob(i,lfedge(j,iedge)),((jj-1)*3)+k)=sol(jj,k)*safety;
                                end
                            end
                        end
                    end %ed
                end
            end % fac
        end % elem
        
        
        % loop over elemets in mesh
        for i = 1:nelem %elem
            % local face to edge connections
            % lfedge contains the edge numbers on each face
            if eltype(i) ==1
                lfedge(1,1)=1;
                lfedge(1,2)=5;
                lfedge(1,3)=6;
                lfedge(4,1)=3;
                lfedge(4,2)=2;
                lfedge(4,3)=1;
            else
                lfedge(1,1)=1;
                lfedge(1,2)=6;
                lfedge(1,3)=5;
                lfedge(4,1)=2;
                lfedge(4,2)=3;
                lfedge(4,3)=1;
            end
            lfedge(2,1)=2;
            lfedge(2,2)=4;
            lfedge(2,3)=6;
            
            lfedge(3,1)=3;
            lfedge(3,2)=4;
            lfedge(3,3)=5;
            
            % loop over faces in mesh
            for j=1:4 %fac
                if cond(i,j)==bcno && mat(i)==2
                    %set up vertices of the face and move back to sphere surface
                    if j==1
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
                    
                    % transfer coordinates to local array
                    xy(1:4,1:3)=coord(intma(i,1:4),1:3);
                    
                    af = zeros(order*(order-1)/2,order*(order-1)/2);
                    rf = zeros(order*(order-1)/2,3);
                    
                    if gorder >=2
                        % correct the face
                        for pp=1:nipf
                            l(1)=1.d0-intfxi(pp)-intfeta(pp);
                            l(2)=intfxi(pp);
                            l(3)=intfeta(pp);
                            
                            xit = l*p(:,1);
                            etat = l*p(:,2);
                            zetat = l*p(:,3);
                            
                            % compute H1 basis for this face
                            ph=basish1(gesizet,xit,etat,zetat,gorder+1,eltype(i));
                            
                            % initialize correction term
                            gc = ph(1:4,1)'*xy(1:4,1:3);
                            rd=sqrt((gc(1)^2)+(gc(2)^2)+(gc(3)^2));
                            x=gc(1)/rd*rin;
                            y=gc(2)/rd*rin;
                            z=gc(3)/rd*rin;
                            
                            % integrate in to correction term missing jj loop
                            % in L2 projection check
                            for ii=1:3
                                for k=1:order
                                    for jj=1:3
                                        row=4+lfedge(j,jj)+6*(k-1);
                                        gc(ii)=gc(ii)+(ph(row)*...
                                            edgecof(glob(i,lfedge(j,jj)),((k-1)*3)+ii));
                                    end
                                end
                            end
                            
                            rd=sqrt((gc(1)^2)+(gc(2)^2)+(gc(3)^2));
                            x=gc(1)/rd*rin;
                            y=gc(2)/rd*rin;
                            z=gc(3)/rd*rin;
                            
                            for ii=1:order*(order-1)/2
                                row=4+6*order+(j-1)*order*(order-1)/2+ii;
                                for jj=1:order*(order-1)/2
                                    col=4+6*order+(j-1)*order*(order-1)/2+jj;
                                    af(ii,jj)=af(ii,jj)+(intfw(pp)*ph(row)*ph(col));
                                end
                                rf(ii,1)=rf(ii,1)+(ph(row)*intfw(pp)*(x-gc(1)));
                                rf(ii,2)=rf(ii,2)+(ph(row)*intfw(pp)*(y-gc(2)));
                                rf(ii,3)=rf(ii,3)+(ph(row)*intfw(pp)*(z-gc(3)));
                            end
                        end
                        % solve the system for the g coefficents
                        solf = af\rf;
                        for ii=1:order*(order-1)/2
                            for jj=1:3
                                if facecof(globfa(i,j),(3*(ii-1))+jj)~=0
                                    if abs(facecof(globfa(i,j),(3*(ii-1))+jj))> ...
                                            abs(solf(ii,jj)*safety)
                                        facecof(globfa(i,j),(3*(ii-1))+jj)=solf(ii,jj)*safety;
                                    end
                                else
                                    facecof(globfa(i,j),(3*(ii-1))+jj)=solf(ii,jj)*safety;
                                end
                            end
                        end
                    end
                end
            end % fac
        end%elem
        
        % check the surface area and volume of the sphere
        testarea = 0;
        for i =1:nelem %elem
            % transfer coordinates to local array
            xy = coord(intma(i,1:4),1:3);
            
            flag=0;
            for j=1:6
                for pp=1:gorder
                    for k=1:3
                        lec(j,pp,k)=edgecof(glob(i,j),((pp-1)*3)+k);
                        if abs(edgecof(glob(i,j),((pp-1)*3)+k))~=0
                            flag=1;
                        end
                        % if(abs(edgecof(global(i,j),((p-1)*3)+k)).gt.0.00000001)flag=1
                    end
                end
            end
            
            for j=1:4
                for pp=1:gorder*(gorder-1)/2
                    for k=1:3
                        lfc(j,pp,k)=facecof(globfa(i,j),((pp-1)*3)+k);
                        if facecof(globfa(i,j),((pp-1)*3)+k)~=0
                            flag=1;
                        end
                    end
                end
            end
            
            % loop over faces in mesh
            for j=1:4 % fac
                if cond(i,j) == bcno && mat(i) ==2
                    % loop over integration points pn face correct the face
                    for pp=1:nipf
                        % compute integration point locations
                        l(1)=1.d0-intfxi(pp)-intfeta(pp);
                        l(2)=intfxi(pp);
                        l(3)=intfeta(pp);
                        
                        xit = l*p(1:3,1);
                        etat = l*p(1:3,2);
                        zetat = l*p(1:3,3);
                        
                        area=surdet(eltype(i),xit,etat,zetat,j,xy,lec,lfc,gorder+1);
                        
                        testarea = testarea + abs(area*intfw(pp));
                    end
                end
            end % fac
        end % elem
        
        area = 4*pi*rin^2;
        display(['The relative error in the surface of the sphere is', num2str(abs(area-testarea)/area)])
        
        % check also volume
        testarea =0;
        for i=1:nelem % elem
               xy = coord(intma(i,1:4),1:3);
                
                flag=0;
                for j=1:6
                    for pp=1:gorder
                        for k=1:3
                            lec(j,pp,k,i)=edgecof(glob(i,j),((pp-1)*3)+k);
                            if abs(edgecof(glob(i,j),((pp-1)*3)+k))~=0
                                flag=1;
                            end
                            % if(abs(edgecof(global(i,j),((p-1)*3)+k)).gt.0.00000001)flag=1
                        end
                    end
                end
                
                for j=1:4
                    for pp=1:gorder*(gorder-1)/2
                        for k=1:3
                            lfc(j,pp,k,i)=facecof(globfa(i,j),((pp-1)*3)+k);
                            if facecof(globfa(i,j),((pp-1)*3)+k)~=0
                                flag=1;
                            end
                        end
                    end
                end
                
                
                
                gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
                
                mycoord = zeros(gesizet,3);
                
                %transfer coefficents to locations vertices
                mycoord(1:4,1:3,i) = xy(1:4,1:3);
                
                % edge functions
                for ii=1:6
                    for p=1:gorder
                        for j=1:3
                            
                            mycoord(4+ii+6*(p-1),j,i)=lec(ii,p,j,i);
                        end
                    end
                end
                
                % face functions
                for iii=1:4
                    for ii=1:(gorder-1)*gorder/2
                        for j=1:3
                            mycoord(4+6*gorder+(iii-1)*gorder*(gorder-1)/2+ii,j,i)= lfc(iii,ii,j,i);
                        end
                    end
                end
            if mat(i) ==2
                % transfer coordinates to local array
             
                
                if eltype(i)==1
                    gphx=gphx1;
                    gphy=gphy1;
                    gphz=gphz1;
                else
                    gphx=gphx2;
                    gphy=gphy2;
                    gphz=gphz2;
                end
                
                for pp = 1:nipv
                    
                    gph(1:gesizet,1)=gphx(pp,1:gesizet)';
                    gph(1:gesizet,2)=gphy(pp,1:gesizet)';
                    gph(1:gesizet,3)=gphz(pp,1:gesizet)';
                    
                    %                 [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian(xy,intxi(pp),inteta(pp),intzeta(pp),...
                    %                                                     eltype(i),lec,lfc,flag,gorder,mycoord);
                    [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord(:,:,i)) ;
                    
                    testarea = testarea + det*intw(pp);
                end
            end
        end % elem
        
        area = 4/3*pi*rin^3;
        display(['The relative error in the volume of the sphere is', num2str(abs(area-testarea)/area)])
        
        if probdata.jb.gag(in+1)~=1 %str2num(jobdata{19+in*3})
            break
        end
        in = in+1;
    end
    
end
mesh.Coordinates=coord;
mesh.mycoord=mycoord;
mesh.lec=lec;
mesh.lfc=lfc;
