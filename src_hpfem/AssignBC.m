 % Function creates a unique numbering of the faces and builds the element
% to face number connectivity.

function [mesh]=AssignBC(mesh,ProblemData)

% iface = the node numbers, element number and boundary data
% globfa = contains the local to global face connectivities
% cond = the boundary type for faces
% bcedge = the boundary type for edges

% Extract relevant data form mesh structure
npoin=mesh.npoin;
nelem=mesh.Nelements;
nboun=mesh.nboun;
nside=mesh.edge.nside;
ielem=mesh.intma;
iboun=mesh.bsido;
glob=mesh.edge.glob;
bcosurf=mesh.bcosurf;


% Define boundary conditions for faces based on face labelling from mesh
bcosurf(1)=ProblemData.bctypeouter;    % Dirichlet or Neumann as defined in problem file
bcosurf(2)=5;                          % Coil interface 
bcosurf(3)=4;                          % Conductor interface
bcosurf(4)=6;                          % Conductor interface (Dirichlet for mechanics)
bcosurf(7)=7;                          % Symmetry BC (x plane)
bcosurf(8)=8;                          % Symmetry BC (y plane)
bcosurf(9)=3;                          % Neumann BC
bcosurf(10)=10;                        % Symmetry EM and antisymmetry Mechanics BC (x plane)
bcosurf(12)=12;                        % Antisymmetry EM and symmetry Mechanics BC  (z plane)
bcosurf(13)=13;                        % Symmetry EM and antisymmetry mechanics (z plane)

% Save this to mesh structure
mesh.bcosurf=bcosurf;

nfaceold = 2*nelem+nboun/2;
no = [1, 1, 1, 4 
      2, 4, 3, 3
      3, 2, 4, 2
      4, 3, 2, 1];
% no = [4, 1, 1, 1  
%       3, 3, 4, 2
%       2, 4, 2, 3 ];
  
iface = zeros(nfaceold,5);
globfa = zeros(nelem,4);
markf = zeros(1,2*npoin);
mark1 = zeros(1,2*npoin);
er = zeros(1,16);

% assume that all surfaces are type 1 bcondition
% bsido(i,5) = global surface number (6 surfaces for cube)
% set up boundcond
cond = zeros(nelem,4);
bcedge = zeros(1,nside);
bcedge2 = zeros(1,nside);
bcedge3 = zeros(1,nside);
bcedge4 = zeros(1,nside);
bcedge5 = zeros(1,nside);
bcedge6 = zeros(1,nside);
bcedge7 = zeros(1,nside);
bcedge8 = zeros(1,nside);
bcvertex = zeros(1,npoin);
bcvertex2 = zeros(1,npoin);
bcvertex3 = zeros(1,npoin);
bcvertex4 = zeros(1,npoin);
bcvertex5 = zeros(1,npoin);
bcvertex6 = zeros(1,npoin);
bcvertex7 = zeros(1,npoin);
bcvertex8 = zeros(1,npoin);
iel = 0; 

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6)); 
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=6 
                 
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
                  
               end
               if bcedge(glob(iele,2)) ~=6 
                   
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
             
               end
               if bcedge(glob(iele,3)) ~=6 
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                  
               end
            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,3)) ~=6 
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,4)) ~=6 
                  
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,5)) ~=6 
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
                 
               end
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,2)) ~=6 
                  
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
                 
               end
               if bcedge(glob(iele,4)) ~=6
                   
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,6)) ~=6
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=6
                   
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
               end
               if bcedge(glob(iele,5)) ~=6
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
               
               end
               if bcedge(glob(iele,6)) ~=6
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end


% Create secondary edge number to find edges shared between two faces
for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                 
                if cond(iele,4)==12
                    bcedge2(glob(iele,1))=8;
                  
                   
                    bcedge2(glob(iele,2))=8;
             
                   
                    bcedge2(glob(iele,3))=8;
                end
                  

            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                   
                if cond(iele,3)==12
                    bcedge2(glob(iele,3))=8;
                   
                  
                    bcedge2(glob(iele,4))=8;

                  
                    bcedge2(glob(iele,5))=8;
                end
                 
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,2)==12
                  
                    bcedge2(glob(iele,2))=8;
                 
                   
                    bcedge2(glob(iele,4))=8;
                   
                   
                    bcedge2(glob(iele,6))=8;
                end
                   
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,1)==12
                   
                    bcedge2(glob(iele,1))=8;
                  
                    bcedge2(glob(iele,5))=8;
                   
                    bcedge2(glob(iele,6))=8;
                end
                   
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end



% Third edge numbering for finding shared edges
for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                 
                if cond(iele,4)==8
                    bcedge3(glob(iele,1))=6;
                  
                   
                    bcedge3(glob(iele,2))=6;
             
                   
                    bcedge3(glob(iele,3))=6;
                end
                  

            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                   
                if cond(iele,3)==8
                    bcedge3(glob(iele,3))=6;
                   
                  
                    bcedge3(glob(iele,4))=6;

                  
                    bcedge3(glob(iele,5))=6;
                end
                 
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,2)==8
                  
                    bcedge3(glob(iele,2))=6;
                 
                   
                    bcedge3(glob(iele,4))=6;
                   
                   
                    bcedge3(glob(iele,6))=6;
                end
                   
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,1)==8
                   
                    bcedge3(glob(iele,1))=6;
                  
                    bcedge3(glob(iele,5))=6;
                   
                    bcedge3(glob(iele,6))=6;
                end
                   
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

% Third edge numbering for finding shared edges
for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                 
                if cond(iele,4)==7
                    bcedge8(glob(iele,1))=6;
                  
                   
                    bcedge8(glob(iele,2))=6;
             
                   
                    bcedge8(glob(iele,3))=6;
                end
                  

            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                   
                if cond(iele,3)==7
                    bcedge8(glob(iele,3))=6;
                   
                  
                    bcedge8(glob(iele,4))=6;

                  
                    bcedge8(glob(iele,5))=6;
                end
                 
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,2)==7
                  
                    bcedge8(glob(iele,2))=6;
                 
                   
                    bcedge8(glob(iele,4))=6;
                   
                   
                    bcedge8(glob(iele,6))=6;
                end
                   
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,1)==7
                   
                    bcedge8(glob(iele,1))=6;
                  
                    bcedge8(glob(iele,5))=6;
                   
                    bcedge8(glob(iele,6))=6;
                end
                   
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                 
                if cond(iele,4)==10
                    bcedge4(glob(iele,1))=6;
                  
                   
                    bcedge4(glob(iele,2))=6;
             
                   
                    bcedge4(glob(iele,3))=6;
                end
                  

            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                   
                if cond(iele,3)==10
                    bcedge4(glob(iele,3))=6;
                   
                  
                    bcedge4(glob(iele,4))=6;

                  
                    bcedge4(glob(iele,5))=6;
                end
                 
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,2)==10
                  
                    bcedge4(glob(iele,2))=6;
                 
                   
                    bcedge4(glob(iele,4))=6;
                   
                   
                    bcedge4(glob(iele,6))=6;
                end
                   
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,1)==10
                   
                    bcedge4(glob(iele,1))=6;
                  
                    bcedge4(glob(iele,5))=6;
                   
                    bcedge4(glob(iele,6))=6;
                end
                   
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

% Fourth edge numbering for finding shared edges
for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                 
                if cond(iele,4)==7
                    bcedge5(glob(iele,1))=13;
                  
                   
                    bcedge5(glob(iele,2))=13;
             
                   
                    bcedge5(glob(iele,3))=13;
                end
                  

            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                   
                if cond(iele,3)==7
                    bcedge5(glob(iele,3))=13;
                   
                  
                    bcedge5(glob(iele,4))=13;

                  
                    bcedge5(glob(iele,5))=13;
                end
                 
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,2)==7
                  
                    bcedge5(glob(iele,2))=13;
                 
                   
                    bcedge5(glob(iele,4))=13;
                   
                   
                    bcedge5(glob(iele,6))=13;
                end
                   
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,1)==7
                   
                    bcedge5(glob(iele,1))=13;
                  
                    bcedge5(glob(iele,5))=13;
                   
                    bcedge5(glob(iele,6))=13;
                end
                   
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                 
                if cond(iele,4)==8
                    bcedge6(glob(iele,1))=13;
                  
                   
                    bcedge6(glob(iele,2))=13;
             
                   
                    bcedge6(glob(iele,3))=13;
                end
                  

            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                   
                if cond(iele,3)==8
                    bcedge6(glob(iele,3))=13;
                   
                  
                    bcedge6(glob(iele,4))=13;

                  
                    bcedge6(glob(iele,5))=13;
                end
                 
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,2)==8
                  
                    bcedge6(glob(iele,2))=13;
                 
                   
                    bcedge6(glob(iele,4))=13;
                   
                   
                    bcedge6(glob(iele,6))=13;
                end
                   
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,1)==8
                   
                    bcedge6(glob(iele,1))=13;
                  
                    bcedge6(glob(iele,5))=13;
                   
                    bcedge6(glob(iele,6))=13;
                end
                   
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                 
                if cond(iele,4)==12
                    bcedge7(glob(iele,1))=6;
                  
                   
                    bcedge7(glob(iele,2))=6;
             
                   
                    bcedge7(glob(iele,3))=6;
                end
                  

            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                   
                if cond(iele,3)==12
                    bcedge7(glob(iele,3))=6;
                   
                  
                    bcedge7(glob(iele,4))=6;

                  
                    bcedge7(glob(iele,5))=6;
                end
                 
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,2)==12
                  
                    bcedge7(glob(iele,2))=6;
                 
                   
                    bcedge7(glob(iele,4))=6;
                   
                   
                    bcedge7(glob(iele,6))=6;
                end
                   
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
                if cond(iele,1)==12
                   
                    bcedge7(glob(iele,1))=6;
                  
                    bcedge7(glob(iele,5))=6;
                   
                    bcedge7(glob(iele,6))=6;
                end
                   
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6
                 
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
                  
               end
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6
                   
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
             
               end
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                  
               end
            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6
                  
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
                 
               end
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6
                  
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
                 
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6
                   
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6
                   
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
               
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end


for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=6
                 
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
                  
               end
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=6
                   
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
             
               end
               if  bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=6
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                  
               end
            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=6
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=6
                  
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if  bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=6
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
                 
               end
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=6
                  
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
                 
               end
               if  bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=6
                   
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if  bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=6
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=6
                   
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
               end
               if  bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=6
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
               
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=6
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=8 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=7
                 
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
                  
               end
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=8 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=7
                   
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
             
               end
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=8 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=7
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                  
               end
            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=8 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=7
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                   
               end
               if  bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=8 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=7
                  
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if  bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=8 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=7
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
                 
               end
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=8 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=7
                  
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
                 
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=8 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=7
                   
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=8 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=7
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=8 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=7
                   
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
               end
               if  bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=8 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=7
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
               
               end
               if  bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=8 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=7
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=10 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=8
                 
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
                  
               end
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=10 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=8
                   
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
             
               end
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=10 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=8
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                  
               end
            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=10 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=8
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=10 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=8
                  
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=10 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=8
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
                 
               end
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=10 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=8
                  
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
                 
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=10 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=8
                   
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=10 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=8
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=10 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=8
                   
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=10 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=8
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
               
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=10 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=8
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end


for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if  bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=12 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=8 && bcedge(glob(iele,1)) ~=10
                 
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
                  
               end
               if  bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=12 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=8 && bcedge(glob(iele,2)) ~=10
                   
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
             
               end
               if  bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=12 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=8 && bcedge(glob(iele,3)) ~=10
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                  
               end
            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=12 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=8 && bcedge(glob(iele,3)) ~=10
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=12 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=8 && bcedge(glob(iele,4)) ~=10
                  
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=12 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=8 && bcedge(glob(iele,5)) ~=10
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
                 
               end
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=12 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=8 && bcedge(glob(iele,2)) ~=10
                  
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
                 
               end
               if  bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=12 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=8 && bcedge(glob(iele,4)) ~=10
                   
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=12 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=8 && bcedge(glob(iele,6)) ~=10
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=12 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=8 && bcedge(glob(iele,1)) ~=10
                   
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=12 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=8 && bcedge(glob(iele,5)) ~=10
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
               
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=12 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=8 && bcedge(glob(iele,6)) ~=10
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    iele = iboun(i,4);
    is = 0;
    for j=1:4
        io1 = ielem(iele,no(1,j));
        io2 = ielem(iele,no(2,j));
        io3 = ielem(iele,no(3,j));
        if ((io1==ip1 && io2==ip2 && io3==ip3) || (io1==ip2 && io2==ip3 && io3==ip1) ...
            || (io1==ip3 && io2==ip1 && io3==ip2) || (io1==ip1 && io2==ip3 && io3==ip2) ...
            || (io1==ip2 && io2==ip1 && io3==ip3) || (io1==ip3 && io2==ip2 && io3==ip1))
            iel = iel+1;
            is = is+1;
            if j==1
                cond(iele,4) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=2 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=8 && bcedge(glob(iele,1)) ~=10 && bcedge(glob(iele,1)) ~=12
                 
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
                  
               end
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=2 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=8 && bcedge(glob(iele,2)) ~=10 && bcedge(glob(iele,2)) ~=12
                   
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
             
               end
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=2 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=8 && bcedge(glob(iele,3)) ~=10 && bcedge(glob(iele,3)) ~=12
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                  
               end
            elseif j==2
                cond(iele,3) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,3)) ~=13 && bcedge(glob(iele,3)) ~=6 && bcedge(glob(iele,3)) ~=2 && bcedge(glob(iele,3)) ~=7 && bcedge(glob(iele,3)) ~=8 && bcedge(glob(iele,3)) ~=10 && bcedge(glob(iele,3)) ~=12
                   
                    bcedge(glob(iele,3))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=2 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=8 && bcedge(glob(iele,4)) ~=10 && bcedge(glob(iele,4)) ~=12
                  
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=2 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=8 && bcedge(glob(iele,5)) ~=10 && bcedge(glob(iele,5)) ~=12
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
                 
               end
            elseif j==3
                cond(iele,2) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,2)) ~=13 && bcedge(glob(iele,2)) ~=6 && bcedge(glob(iele,2)) ~=2 && bcedge(glob(iele,2)) ~=7 && bcedge(glob(iele,2)) ~=8 && bcedge(glob(iele,2)) ~=10 && bcedge(glob(iele,2)) ~=12
                  
                    bcedge(glob(iele,2))=bcosurf(iboun(i,6));
                 
               end
               if bcedge(glob(iele,4)) ~=13 && bcedge(glob(iele,4)) ~=6 && bcedge(glob(iele,4)) ~=2 && bcedge(glob(iele,4)) ~=7 && bcedge(glob(iele,4)) ~=8 && bcedge(glob(iele,4)) ~=10 && bcedge(glob(iele,4)) ~=12
                   
                    bcedge(glob(iele,4))=bcosurf(iboun(i,6));
                   
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=2 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=8 && bcedge(glob(iele,6)) ~=10 && bcedge(glob(iele,6)) ~=12
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            else
                cond(iele,1) = bcosurf(iboun(i,6));
                %boundary conditions on edges
               if bcedge(glob(iele,1)) ~=13 && bcedge(glob(iele,1)) ~=6 && bcedge(glob(iele,1)) ~=2 && bcedge(glob(iele,1)) ~=7 && bcedge(glob(iele,1)) ~=8 && bcedge(glob(iele,1)) ~=10 && bcedge(glob(iele,1)) ~=12
                   
                    bcedge(glob(iele,1))=bcosurf(iboun(i,6));
               end
               if bcedge(glob(iele,5)) ~=13 && bcedge(glob(iele,5)) ~=6 && bcedge(glob(iele,5)) ~=2 && bcedge(glob(iele,5)) ~=7 && bcedge(glob(iele,5)) ~=8 && bcedge(glob(iele,5)) ~=10 && bcedge(glob(iele,5)) ~=12
                  
                  
                    bcedge(glob(iele,5))=bcosurf(iboun(i,6));
               
               end
               if bcedge(glob(iele,6)) ~=13 && bcedge(glob(iele,6)) ~=6 && bcedge(glob(iele,6)) ~=2 && bcedge(glob(iele,6)) ~=7 && bcedge(glob(iele,6)) ~=8 && bcedge(glob(iele,6)) ~=10 && bcedge(glob(iele,6)) ~=12
                   
                    bcedge(glob(iele,6))=bcosurf(iboun(i,6));
                   
               end
            end
        end
    end
    if is==0 || is >1
        display('no bc face found for element')
    end
end



% Find the edges shared between conflicting faces and re-number them
for ii=1:length(bcedge)
    if bcedge(ii)==bcedge2(ii) && bcedge(ii) ~= 0
        bcedge(ii)=128;
    end
    if bcedge(ii)==bcedge3(ii) && bcedge(ii) ~= 0
        bcedge(ii)=68;
    end
    if bcedge(ii)==bcedge4(ii) && bcedge(ii) ~= 0
        bcedge(ii)=610;
    end
    if bcedge(ii)==bcedge5(ii) && bcedge(ii) ~= 0
        bcedge(ii)=137;
    end
    if bcedge(ii)==bcedge6(ii) && bcedge(ii) ~= 0
        bcedge(ii)=138;
    end 
    if bcedge(ii)==bcedge7(ii) && bcedge(ii) ~= 0
        bcedge(ii)=612;
    end 
    if bcedge(ii)==bcedge8(ii) && bcedge(ii) ~= 0
        bcedge(ii)=67;
    end 
end

ihelp(1:npoin) = zeros(1,npoin);
for ib=1:nboun  %20
    ihelp(iboun(ib,1)) =1;
    ihelp(iboun(ib,2)) =1;
    ihelp(iboun(ib,3)) =1;
end %20
for i = 1:nfaceold %40
    iface(i,4) =0;
    iface(i,5) =0;
end %40

% find number of elements surrounding a node
 markf(1:npoin) = zeros(1,npoin);
for in = 1:4  %70
    for ie=1:nelem % 60
        ip=ielem(ie,in);
        markf(ip) = markf(ip)+1;
    end %60
end %70

for i=1:npoin %75
    if markf(ip)==0
        er(1) = er(1)+1;
        display('point does not appear in the connectivities')
        continue
    end
    if ihelp(ip) ==1 
        continue
    end
    if mod(markf(ip),2) ~=0
        er(2) = er(2)+1;
        display('wrong number of element arround point')
    end
end %75

mark1(1) = 0;
for ip = 2:npoin % 80
    mark1(ip) = mark1(ip-1)+markf(ip-1);
end % 80

markf(1:npoin) = zeros(1,npoin);

% find the elements surrounding each node  
for in=1:4 % 110
    for ie =1:nelem %100
        ip = ielem(ie,in);
        markf(ip) = markf(ip)+1;
        jloca = mark1(ip)+markf(ip);
        ihelp(jloca) = ie;
    end %100
end %110

% split into faces
iloca =0;
for ip=1:npoin   %240
    iloc1 = iloca;
    iele = markf(ip);
    if iele ==0
        continue;
    end
    iwher = mark1(ip);
    ip1 = ip;
    
    for iel = 1:iele   % 230
        ie = ihelp(iwher+iel);
        for in = 1:4 %120
            in1 = in;
            ipt = ielem(ie,in);
            if ipt == ip
                break;
            end
        end %120
        
        mk1 = 0;
        mk2 = 0;
        mk3 = 0;
        for j=1:3  %220
            in2 = in1+j;
            if in2 >4
                in2 = in2-4;
            end
            if in2 == in1
                in2 = in2+1;
            end
            if in2 >4
                in2 = in2-4;
            end
            ip2 = ielem(ie,in2);
            if ip2 <ip1
                continue;
            end
            in3 = in2+1;
            if in3 >4
                in3 = in3 -4;
            end
            if in3 ==in1
                in3 = in3 +1;
            end
            if in3 >4
                in3 = in3 -4;
            end
            ip3 = ielem(ie,in3);
            if ip3 <ip1
                continue;
            end
            if iloca == iloc1 
                mk1 = 150;
            end
            
            if mk1 ~=150
            for is = (iloc1+1):iloca % 140
                jloca = is;
                if (iface(is,2)==ip2 && iface(is,3)==ip3) || (iface(is,2)==ip3 && iface(is,3)==ip2)
                    mk2 = 180;
                    break
                end
                mk3 = is;
            end %140
            end
            
            mk4 = 0;
            mk5 = 0;
            if mk1 == 150 || mk3==is
                % new face
                iloca = iloca +1;
                for ko = 1:4 %160
                    io1=ielem(ie,no(1,ko));
                    io2=ielem(ie,no(2,ko));
                    io3=ielem(ie,no(3,ko));
                    if (ip1==io1 || ip1==io2 || ip1==io3) && ...
                       (ip2==io1 || ip2==io2 || ip2==io3) && ...
                       (ip3==io1 || ip3==io2 || ip3==io3)
                        if ip1 ==io1
                            iface(iloca,1) = io1;
                            iface(iloca,2) = io2;
                            iface(iloca,3) = io3;
                        elseif ip1 == io2
                            iface(iloca,1) = io2;
                            iface(iloca,2) = io3;
                            iface(iloca,3) = io1;
                        else
                            iface(iloca,1) = io3;
                            iface(iloca,2) = io1;
                            iface(iloca,3) = io2;
                        end
                        % locate inside globfa using Joes numbering
                        if ko==1
                            globfa(ie,4) = iloca;
                        elseif ko==2
                            globfa(ie,3) = iloca;
                        elseif ko==3
                            globfa(ie,2) = iloca;
                        else
                            globfa(ie,1) = iloca;
                        end
                        ipos = 4;
                        break
                    end
                    mk4= mk4+1;
                end %160
                if mk4 ==4
                    error('error');
                end
                if iface(iloca,ipos) ~=0
                    error(message('possible dimension error'));
                else
                    iface(iloca,ipos) =ie;
                end

            elseif mk2 == 180
                %old face
                for ko = 1:4 %190
                    io1 = ielem(ie,no(1,ko));
                    io2 = ielem(ie,no(2,ko));
                    io3 = ielem(ie,no(3,ko));
                    if (ip1==io1||ip1==io2||ip1==io3)&&...
                            (ip2==io1||ip2==io2||ip2==io3)&&...
                            (ip3==io1||ip3==io2||ip3==io3)
                        %store face data in Joe`s Format
                        if ko==1
                            globfa(ie,4) = jloca;
                        elseif ko==2
                            globfa(ie,3) = jloca;
                        elseif ko==3
                            globfa(ie,2) = jloca;
                        else
                            globfa(ie,1) = jloca;
                        end
                        
                        % suspend warnings due to Joe`s new`s orientation
                        if ip1==io1
                            if iface(jloca,1)~=io1||iface(jloca,2)~=io3||iface(jloca,3)~=io2
                                er(16) = er(16)+1;
                            end
                        elseif ip1==io2
                            if iface(jloca,1)~=io2||iface(jloca,2)~=io1||iface(jloca,3)~=io3
                                er(16) = er(16)+1;
                            end
                        elseif ip1==io3
                            if iface(jloca,1)~=io3||iface(jloca,2)~=io2||iface(jloca,3)~=io1
                                er(16) = er(16)+1;
                            end
                        else
                            er(3)=er(3)+1;
                            error(message('an old face is badly orriantated'));
                        end
                        ipos = 5;
                        break
                    end
                    mk5 = mk5+1;
                end %190
                if mk5 ==5
                    er(4)=er(4)+1;
                    error(message('an old face is not found'));
                end
                if iface(jloca,ipos) ~=0
                    display('possible dimension error')
%                     error(message('possible dimension error'));
                else
                    iface(jloca,ipos) = ie;
                end
            end    
           
        end % 220  
    end %230
end  %240

if iloca ~= nfaceold
    display('nface ~= iloca due to internal interface')
%    nfaceold,iloca
end
nface = iloca;

%--------------------------------------------------------------------------
for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
     if bcvertex(ip1)~=6
       
    bcvertex(ip1)= bcosurf(iboun(i,6));
       
     end
    if bcvertex(ip2)~=6
        
    bcvertex(ip2)= bcosurf(iboun(i,6));
     
     end
     if bcvertex(ip3)~=6
   
    bcvertex(ip3)= bcosurf(iboun(i,6));
        
     end
end

% Create second, third and fourth vertex numbering to find shared edges


for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcosurf(iboun(i,6))==12
       
    bcvertex2(ip1)= 8;
               
    bcvertex2(ip2)= 8;
  
    bcvertex2(ip3)= 8;
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcosurf(iboun(i,6))==8
       
    bcvertex3(ip1)= 6;
               
    bcvertex3(ip2)= 6;
  
    bcvertex3(ip3)= 6;
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcosurf(iboun(i,6))==7
       
    bcvertex8(ip1)= 6;
               
    bcvertex8(ip2)= 6;
  
    bcvertex8(ip3)= 6;
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcosurf(iboun(i,6))==10
       
    bcvertex4(ip1)= 6;
               
    bcvertex4(ip2)= 6;
  
    bcvertex4(ip3)= 6;
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcosurf(iboun(i,6))==7
       
    bcvertex5(ip1)= 13;
               
    bcvertex5(ip2)= 13;
  
    bcvertex5(ip3)= 13;
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcosurf(iboun(i,6))==8
       
    bcvertex6(ip1)= 13;
               
    bcvertex6(ip2)= 13;
  
    bcvertex6(ip3)= 13;
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcosurf(iboun(i,6))==12
       
    bcvertex7(ip1)= 6;
               
    bcvertex7(ip2)= 6;
  
    bcvertex7(ip3)= 6;
        
    end
end
for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if  bcvertex(ip1)~=13 && bcvertex(ip1)~=8 && bcvertex(ip1)~=7 && bcvertex(ip1)~=6
       
    bcvertex(ip1)= bcosurf(iboun(i,6));
       
    end
    if bcvertex(ip1)~=13 && bcvertex(ip2)~=8 && bcvertex(ip2)~=7 && bcvertex(ip2)~=6
        
    bcvertex(ip2)= bcosurf(iboun(i,6));
     
    end
    if bcvertex(ip3)~=13 && bcvertex(ip3)~=8 && bcvertex(ip3)~=7 && bcvertex(ip3)~=6
  
    bcvertex(ip3)= bcosurf(iboun(i,6));
        
    end
end
for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcvertex(ip1)~=13 && bcvertex(ip1)~=6 && bcvertex(ip1)~=7 && bcvertex(ip1)~=8 && bcvertex(ip1)~=10
       
    bcvertex(ip1)= bcosurf(iboun(i,6));
       
    end
    if bcvertex(ip2)~=13 && bcvertex(ip2)~=6 && bcvertex(ip2)~=7 && bcvertex(ip2)~=8 && bcvertex(ip2)~=10
        
    bcvertex(ip2)= bcosurf(iboun(i,6));
     
    end
    if bcvertex(ip3)~=13 && bcvertex(ip3)~=6 && bcvertex(ip3)~=7 && bcvertex(ip3)~=8 && bcvertex(ip3)~=10
  
    bcvertex(ip3)= bcosurf(iboun(i,6));
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcvertex(ip1)~=13 && bcvertex(ip1)~=6 && bcvertex(ip1)~=7 && bcvertex(ip1)~=8 && bcvertex(ip1)~=10 && bcvertex(ip1)~=12
       
    bcvertex(ip1)= bcosurf(iboun(i,6));
       
    end
    if bcvertex(ip2)~=13 && bcvertex(ip2)~=6 && bcvertex(ip2)~=7 && bcvertex(ip2)~=8 && bcvertex(ip2)~=10 && bcvertex(ip2)~=12
        
    bcvertex(ip2)= bcosurf(iboun(i,6));
     
    end
    if bcvertex(ip3)~=13 && bcvertex(ip3)~=6 && bcvertex(ip3)~=7 && bcvertex(ip3)~=8 && bcvertex(ip3)~=10 && bcvertex(ip3)~=12
  
    bcvertex(ip3)= bcosurf(iboun(i,6));
        
    end
end

for i = 1:nboun
    ip1 = iboun(i,1);
    ip2 = iboun(i,2);
    ip3 = iboun(i,3);
    if bcvertex(ip1)~=13 && bcvertex(ip1)~=6 && bcvertex(ip1)~=7 && bcvertex(ip1)~=8 && bcvertex(ip1)~=2 && bcvertex(ip1)~=10 && bcvertex(ip1)~=12
       
    bcvertex(ip1)= bcosurf(iboun(i,6));
       
    end
    if bcvertex(ip2)~=13 && bcvertex(ip2)~=6 && bcvertex(ip2)~=7 && bcvertex(ip2)~=8 && bcvertex(ip2)~=2 && bcvertex(ip2)~=10 && bcvertex(ip2)~=12
        
    bcvertex(ip2)= bcosurf(iboun(i,6));
     
    end
    if bcvertex(ip3)~=13 && bcvertex(ip3)~=6 && bcvertex(ip3)~=7 && bcvertex(ip3)~=8 && bcvertex(ip3)~=2 && bcvertex(ip3)~=10 && bcvertex(ip3)~=12
  
    bcvertex(ip3)= bcosurf(iboun(i,6));
        
    end
end

% Find the shared edges and re-number them
for ii=1:length(bcvertex)
    if bcvertex(ii)==bcvertex2(ii) && bcvertex(ii) ~= 0
        bcvertex(ii)=128;
    end
    if bcvertex(ii)==bcvertex3(ii) && bcvertex(ii) ~= 0
        bcvertex(ii)=68;
    end
    if bcvertex(ii)==bcvertex4(ii) && bcvertex(ii) ~= 0
        bcvertex(ii)=610;
    end
    if bcvertex(ii)==bcvertex5(ii) && bcvertex(ii) ~= 0
        bcvertex(ii)=137;
    end
    if bcvertex(ii)==bcvertex6(ii) && bcvertex(ii) ~= 0
        bcvertex(ii)=138;
    end
    if bcvertex(ii)==bcvertex7(ii) && bcvertex(ii) ~= 0
        bcvertex(ii)=612;
    end
    if bcvertex(ii)==bcvertex8(ii) && bcvertex(ii) ~= 0
        bcvertex(ii)=67;
    end
end


mesh.face.globfa=globfa;
mesh.face.cond=cond;
mesh.face.nface=nface;
mesh.edge.bcedge=bcedge;
mesh.bcvertex=bcvertex;
    
