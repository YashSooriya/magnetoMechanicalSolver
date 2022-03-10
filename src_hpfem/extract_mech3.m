function [Xzz,XEgEg,XFgFg, XFF, XIgIg, XII, nz, nheg, nhfg, nhf,nhig, nhi,stiffVV,stiffEE,stiffFF,stiffII,nVV,nEE,nFF,nII]=extract_mech3(unknown,Kpre,K,ProblemData)

%=========================================================================
% Extract data from structures
%=========================================================================
order= ProblemData.order;
orderH1=ProblemData.orderH1;
nunkVertex=unknown.Mech.nunkVertex;
nunkEdges=unknown.Mech.nunkEdges;
nunkFaces=unknown.Mech.nunkFaces;
nunkInteriors=unknown.Mech.nunkInteriors;
unkz=unknown.EM.unkz;
unksid=unknown.EM.unksid;
unkfatp1=unknown.EM.unkfatp1;
unkfatp2=unknown.EM.unkfatp2;
unkfatp3=unknown.EM.unkfatp3;
unkint=unknown.EM.unkint;
unkEdgesMech=unknown.system.unkedgesSystem;
unkFacesMech=unknown.system.unkfacesSystem;
unkInteriorsMech=unknown.system.unkinteriorsSystem;

%-------------------------------------------------------------------------
% Initialise matrices
%-------------------------------------------------------------------------
stiffEE=[];
stiffFF=[];
stiffII=[];
nEE=0;
nFF=0;
nII=0;

%--------------------------------------------------------------------------
% Extract the different blocks
%--------------------------------------------------------------------------
% Zero block

% Find number of zero unknowns
nstart=1;
nend=max(unkz);
nz=nend;

% extract zero block
Xzz=K(nstart:nend,nstart:nend);
size(Xzz);

% Find number of edge unknowns
XEgEg=[];
XFgFg=[];
XFF=[];
XIgIg=[];
XII=[];
nheg=0;
nhfg=0;
nhf=0;
nhig=0;
nhi=0;
if order > 0
    
    nstart=nend+1;
    nend=max(max(unksid));
    if nend> 0
        nheg=nend-nz;
        
        % extract zero block
        XEgEg=Kpre(nstart:nend,nstart:nend);

        [XEgEg]=extractsparse(XEgEg,unksid,nstart-1);
    else
        nend=nstart-1;
        % no edge gradients!
    end
end

if order > 1
    % Find number of gradient face unknowns
    nstart=nend+1;
    nend=max(max(unkfatp1));
    if nend > 0
        nhfg=nend-nz-nheg;
        
        % extract face gradient block
        XFgFg=Kpre(nstart:nend,nstart:nend);
        [XFgFg]=extractsparse(XFgFg,unkfatp1,nstart-1);
        
    else
        nend=nstart-1;
        % no face gradients
    end
    
    % Find number of face (non-gradient) unknowns
    nstart=nend+1;
    nend=max(max(unkfatp3));
    nhf=nend-nz-nheg-nhfg;
    
    % extract non-gradient block
    XFF=Kpre(nstart:nend,nstart:nend);
    % we need to put unkfatp2 and unkfatp3 together!
    unkfatpng=[unkfatp2 unkfatp3];
    
    [XFF]=extractsparse(XFF,unkfatpng,nstart-1);
end

if order > 2
    k=0;
    for ii=0:order-3
        for jj=0:order-3
            for kk=0:order-3
                if ii+jj+kk <= order-3
                    k=k+1;
                end
            end
        end
    end
    nintbas=k;
    
    % Find number of gradient interior unknowns
    nstart=nend+1;
    nend=max(max(unkint(:,1:nintbas)));
    if nend > 0
        nhig=nend-nz-nheg-nhfg-nhf;
        
        % extract gradient interior block
        XIgIg=Kpre(nstart:nend,nstart:nend);
        % block diagonal already
        
    else
        nend=nstart-1;
    end
    
    % Find number of non-gradient interior unknowns
    nstart=nend+1;
    nend=max(max(unkint));
    nhi=nend-nz-nheg-nhfg-nhf-nhig;
    
    % extract interior block
    XII=Kpre(nstart:nend,nstart:nend);
    size(XII);
    
    % block diagonal already
    
end

%==========================================================================
% Mechanical part
%==========================================================================

nstart=nend+1;
nend=nunkVertex;

stiffVV=K(nstart:nend,nstart:nend);
nVV=nend;

if orderH1>1
nstart=nend+1;
nend=nunkEdges;
nEE=nend-nVV;

stiffEE=K(nstart:nend,nstart:nend);
[stiffEE]=extractsparse(stiffEE,unkEdgesMech,nstart-1);
end

if orderH1>2
nstart=nend+1;
nend=nunkFaces;
nFF=nend-nVV-nEE;
stiffFF=K(nstart:nend,nstart:nend);
end
[stiffFF]=extractsparse(stiffFF,unkFacesMech,nstart-1);



if orderH1>3
nstart=nend+1;
nend=nunkInteriors;
stiffII=K(nstart:nend,nstart:nend);
nII=nend-nVV-nFF-nEE;
[stiffII]=extractsparse(stiffII,unkInteriorsMech,nstart-1);
end








