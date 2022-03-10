function [nz,help1,help2,help3]= getnonzero(nelem,nside,esize,mxc,nunkt,unkz,unksid,unkint,order,...
    glob,globfa,unkfatp1,unkfatp2,unkfatp3,helpcoil,maxunk,nef)


% create temporary arrays
nn=6;
temp = zeros(nside,1);
temp2 = zeros(nside,mxc);

for i=1:nelem
    for j=1:nn
        temp(glob(i,j))=temp(glob(i,j))+1;
        if temp(glob(i,j))>mxc
            error(message('increase mxc'));
        end
        temp2(glob(i,j),temp(glob(i,j)))=i;
    end
end

% intialize graph and lp
lp = zeros(nelem,1);
graph = zeros(nelem,mxc);
display('graph structure')

% create graph structure
for i=1:nside %110
    % transfer  elements connected to edge i to local array kp
    if temp(i)>1000
        error(message('increase dimension of kp'));
    end
    for j=1:temp(i) %25
        kp(j)=temp2(i,j);
    end %25
    
    for j=1:temp(i)-1, %100
        nzp=kp(j);
        for k=j+1:temp(i) % 90
            nnp=kp(k);
            for p=1:mxc %60
                % if element already in connectivity list:exit
                coe1=0;
                if graph(nzp,p)==nnp
                    coe1=1;
                    break
                end
                % if some element is in the connectivity list:continue looping
                coe=0;
                if graph(nzp,p)>0
                    coe = p+1;
                    continue
                end
                %  no element found so create this one
                graph(nzp,p)=nnp;
                lp(nzp)=lp(nzp)+1;
                break
            end %60       continue
            if coe1~=1
                if coe ==mxc+1
                    error(message('increase mxc'));
                end
                lp(nnp)=lp(nnp)+1; % 70
                if lp(nnp)<=mxc
                    graph(nnp,lp(nnp))=nzp; %80
                else
                    error(message('increase mxc'));
                end
            end
        end% 90    continue
    end %100  continue
end %110  continue

for i=1:nelem
    lp(i)=lp(i)+1;
    if lp(i)>mxc
        error(message('increase mxc'));
    end
    graph(i,lp(i))=i;
end
 

% compute the number of non-zero enteries for the mass
% and stiffness matrices 
% zero help arrays
help1 = zeros(nef,1);
help2 = zeros(nef,mxc);
      
% help array to flag those enteries already found
help3 = zeros(nef,1);

display('found all linkages')

% obtain non-zero enteries for hcurl basis

% count only face and edge based unknowns
if esize>1000
    error(message('increase out array'));
end

for i=1:nelem
    for j=1:lp(i)
        out=rowfun(unkz,unksid,unkint,order,glob,globfa,graph(i,j),...
                                  unkfatp1,unkfatp2,unkfatp3,esize);

        for p=1:esize % 10
            if out(p)>0
                row=helpcoil(out(p));
            else
                row=0;
            end
            if row<=0 || row>nef
                continue % goto10
            end
            
            % transfer enteries in to help3
            for q=1:help1(row)
                if help2(row,q)==0
                    display('bug');
                end
                help3(help2(row,q))=1;
            end
            
            for q=1:esize % 20
                if out(q)>nunkt
                    display('found bug')
                    %pause
                end
                if out(q)>0
                    col=helpcoil(out(q));
                else
                    col=0;
                end
                
                if col<=0 || col>nef
                    continue %goto20
                end
                
                if help3(col)~=1
                    help1(row)=help1(row)+1;
                    if help1(row)>mxc
                        error(message('increase mxc for help1(row)'));
                    end
                    if row>nef
                        display('found bug2')
                        %pause
                    end
                    help2(row,help1(row))=col;
                    help3(col)=1;
                end
            end %20
            % reintialise help array
            for q=1:help1(row)
                if help2(row,q)>nef
                    display('found bug3')
                    %pause
                end
                help3(help2(row,q))=0;
            end
        end %10
    end
end

% compute actual number of non-zero enteries
nz=0;
for i=1:nef
    for j=1:help1(i)
        % do lower triangle for harwell
        if help2(i,j)>=i
            nz=nz+1;
        end
    end
end

display(['The full mass and stiffness matrices contain',num2str(nef^2),'enteries']);
display(['There are',num2str(nz),'non-zero enteries']);

if nz>maxunk*500
    error(message('nz > = maxunk*500'));
end