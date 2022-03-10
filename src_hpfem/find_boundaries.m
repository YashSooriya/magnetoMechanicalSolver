% Function to automatically detect the internal and external boundaries
% so that user does need to perscribe types to the outer boundary and the 
% object-exterior interface in Netgen.

% Requires object and external boundary not to be connected.

% Any data read in from iboun in mesh file is overwritten.

function [mesh]=find_boundaries(mesh,probdata)

% We first number the faces assume no boundaries have been identified.
% This numbering is then disgarded since we'll do this again once the 
% boundaries are known.

% Extract relevant data from mesh structure
npoin=mesh.npoin;
nelem=mesh.Nelements;
ielem=mesh.intma;
mat=mesh.mat;
subFlag=mesh.subFlag;
matc=mesh.matc;
bsido=mesh.bsido;


%estimate for the number of faces
nfaceold = 4*nelem;
no = [1, 1, 1, 4 
      2, 4, 3, 3
      3, 2, 4, 2
      4, 3, 2, 1];

  
iface = zeros(nfaceold,5);
globfa = zeros(nelem,4);
markf = zeros(1,2*npoin);
markp = zeros(1,2*npoin);
mark1 = zeros(1,2*npoin);
er = zeros(1,16);
ihelp(1:npoin) = zeros(1,npoin);

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
    display('nface ~= iloca due to an internal interface')
%    nfaceold,iloca
end
nface = iloca;

% OK Now we have a list of faces lets identify the boundaries
helpbcface1=zeros(nface,1);
helpbcface2=zeros(nface,2);

% Find which elements share a face
for i=1:nelem
    for j=1:4
        helpbcface1(globfa(i,j))=helpbcface1(globfa(i,j))+1;
        if helpbcface1(globfa(i,j)) <= 2
            % This is a valid number of matching elements and store the
            % element number
            helpbcface2(globfa(i,j),helpbcface1(globfa(i,j)))=i;
        else
            error('too many elements sharing a face')
        end
    end
end

% For each face with only one element we set this an outer boundary

bcosurf(1)=probdata.bctypeouter;%3;%2; % Dirchlet Outer Boundary
           % 3 % Neumann Outer Boundary
bcosurf(2)=5; % Coil interface 
bcosurf(3)=4; % Conductor interface
bcosurf(4)=6;
nboun=0;
iboun=zeros(nface,5); % set this array to be too big since we don't know the number of boundaries
for i=1:nface
    if helpbcface1(i)==1
        % This is an outer boudary
        nboun=nboun+1;
        for j=1:3
            iboun(nboun,j)=iface(i,j); % vertices
        end
        iboun(nboun,4)=helpbcface2(i,1); % element number

            iboun(nboun,5)=1;
            iboun(nboun,6)=bsido(nboun,6);
    end
end

for i=1:nface
    if helpbcface1(i)==2
       % disp('two elements share a face')
       % matc(helpbcface2(i,1)), matc(helpbcface2(i,2))
        if abs(mat(helpbcface2(i,1))-mat(helpbcface2(i,2)) )>=1
            if abs(subFlag(helpbcface2(i,1))-subFlag(helpbcface2(i,2)))==0
        
%        if abs(matc(mat(helpbcface2(i,1)))-matc(mat(helpbcface2(i,2))))==1
            % This is an interface
            nboun=nboun+1;
            for j=1:3
                iboun(nboun,j)=iface(i,j); % vertices
            end
            iboun(nboun,4)=helpbcface2(i,1); % element number
            iboun(nboun,5)=2; % fix 2 as the idendtifier of the interface
            iboun(nboun,6)=bsido(nboun,6);
            nboun=nboun+1;
            for j=1:3
                iboun(nboun,j)=iface(i,j); % vertices
            end
            iboun(nboun,4)=helpbcface2(i,2); % element number
            iboun(nboun,5)=2; % fix 2 as the idendtifier of the interface
            iboun(nboun,6)=bsido(nboun,6);
            %disp('found interface')
           % pause
            end
           if abs(subFlag(helpbcface2(i,1))-subFlag(helpbcface2(i,2)))==2
        
%        if abs(matc(mat(helpbcface2(i,1)))-matc(mat(helpbcface2(i,2))))==1
            % This is an interface
            nboun=nboun+1;
            for j=1:3
                iboun(nboun,j)=iface(i,j); % vertices
            end
            iboun(nboun,4)=helpbcface2(i,1); % element number
            iboun(nboun,5)=3; % fix 2 as the idendtifier of the interface
            iboun(nboun,6)=bsido(nboun,6);
            nboun=nboun+1;
            for j=1:3
                iboun(nboun,j)=iface(i,j); % vertices
            end
            iboun(nboun,4)=helpbcface2(i,2); % element number
            iboun(nboun,5)=3; % fix 2 as the idendtifier of the interface
            iboun(nboun,6)=bsido(nboun,6);
            %disp('found interface')
           % pause
            end
        end
    end
end



mesh.bsido=iboun;
mesh.bcosurf=bcosurf;
mesh.nboun=nboun;

        








    
