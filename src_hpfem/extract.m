% Function to extact out different matrix blocks in order to
% apply the preconditioner.

function [ Xzz,XEgEg,XFgFg, XFF, XIgIg, XII, nz, nheg, nhfg, nhf,nhig, nhi] = extract(unknown,Xpre,order,X )

% Extract relevant data from unknown structure
unkz=unknown.EM.unkz;
unksid=unknown.EM.unksid;
unkfatp1=unknown.EM.unkfatp1;
unkfatp2=unknown.EM.unkfatp2;
unkfatp3=unknown.EM.unkfatp3;
unkint=unknown.EM.unkint;

% Zero block

% Find number of zero unknowns
nstart=1;
nend=max(unkz);
nz=nend;

% extract zero block
%Xzz=Xpre(nstart:nend,nstart:nend);
Xzz=X(nstart:nend,nstart:nend);
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
        XEgEg=Xpre(nstart:nend,nstart:nend);
        %         size(XEgEg)
        %         whos XEgEg
        [XEgEg]=extractsparse(XEgEg,unksid,nstart-1);
        %         size(XEgEg)
        %         whos XEgEg
        %         pause
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
        XFgFg=Xpre(nstart:nend,nstart:nend);
        %         size(XFgFg);
        %         whos XFgFg
        [XFgFg]=extractsparse(XFgFg,unkfatp1,nstart-1);
        %         size(XFgFg)
        %         whos XFgFg
        %         pause
        
    else
        nend=nstart-1;
        % no face gradients
    end
    
    % Find number of face (non-gradient) unknowns
    nstart=nend+1;
    nend=max(max(unkfatp3));
    nhf=nend-nz-nheg-nhfg;
    
    % extract non-gradient block
    XFF=Xpre(nstart:nend,nstart:nend);
    %         size(XFF)
    %         whos XFF
    % we need to put unkfatp2 and unkfatp3 together!
    unkfatpng=[unkfatp2 unkfatp3];
    
    [XFF]=extractsparse(XFF,unkfatpng,nstart-1);
    %         size(XFF)
    %         whos XFF
    %         pause
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
        XIgIg=Xpre(nstart:nend,nstart:nend);
        % block diagonal already
        
    else
        nend=nstart-1;
    end
    
    % Find number of non-gradient interior unknowns
    nstart=nend+1;
    nend=max(max(unkint));
    nhi=nend-nz-nheg-nhfg-nhf-nhig;
    
    % extract interior block
    XII=Xpre(nstart:nend,nstart:nend);
    size(XII);
    
    % block diagonal already
    
end