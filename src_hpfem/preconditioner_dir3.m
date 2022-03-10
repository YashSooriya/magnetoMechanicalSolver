function [out]=preconditioner_dir3(in,arg)
Xzz=arg.Xzz;
nz=arg.nz;

XEgEg=arg.XEgEg;
nheg=arg.nheg;

XFgFg=arg.XFgFg;
nhfg=arg.nhfg;

XFF=arg.XFF;
nhf=arg.nhf;

XIgIg=arg.XIgIg;
nhig=arg.nhig;

XII=arg.XII;
nhi=arg.nhi;
nunkt=arg.nunkt;

% Zero block
inz=in(1:nz);

tmp=arg.L\(arg.P*inz);
tmp2=arg.U\tmp;
outz=arg.Q*tmp2;

%Edge block
if nheg > 0
    inEg=in(nz+1:nz+nheg);
    tmpE=arg.Leg\(arg.Peg*inEg);
    tmpE2=arg.Ueg\tmpE;
    outEg=arg.Qeg*tmpE2;
else
    outEg=[];
end
%Face block
if nhfg > 0
    inFg=in(nz+nheg+1:nz+nheg+nhfg);
    tmpF=arg.Lfg\(arg.Pfg*inFg);
    tmpF2=arg.Ufg\tmpF;
    outFg=arg.Qfg*tmpF2;
else
    outFg=[];
end

%Face block (non-gradients)
if nhf > 0
    inF=in(nz+nheg+nhfg+1:nz+nheg+nhfg+nhf);
    tmpF=arg.Lf\(arg.Pf*inF);
    tmpF2=arg.Uf\tmpF;
    outF=arg.Qf*tmpF2;
else
    outF=[];
end

%Interior block (gradients)
if nhig > 0
    inIg=in(nz+nheg+nhfg+nhf+1:nz+nheg+nhfg+nhf+nhig);
    tmpI=arg.Lig\(arg.Pig*inIg);
    tmpI2=arg.Uig\tmpI;
    outIg=arg.Qig*tmpI2;
else
    outIg=[];
end

%Interior block (non-gradients)
if nhi > 0
    inI=in(nz+nheg+nhfg+nhf+nhig+1:nz+nheg+nhfg+nhf+nhig+nhi);
    tmpI=arg.Li\(arg.Pi*inI);
    tmpI2=arg.Ui\tmpI;
    outI=arg.Qi*tmpI2;
else
    outI=[];
end

%==========================================================================
% Mechanical part
%==========================================================================
XVV=arg.XVV;
nVV=arg.nVV;

XEE=arg.XEE;
nEE=arg.nEE;

XFF=arg.XFF;
nFF=arg.nFF;

XII=arg.XII;
nII=arg.nII;
%If solving just mechanics
%----------------------------------------------------------
% nVV=nVV-nunkt;
%---------------------------------------------------------------


% % Zero block
% inVV=in(1:nVV);
% outVV=pcg(XVV,inVV,1e-8,400);

% Zero block
inVV=in(nunkt+1:nVV);
%inVV=in(1:nVV);

tmp=arg.Lm\(arg.Pm*(inVV));
tmp2=arg.Um\tmp;
outVV=arg.Qm*tmp2;

%outz=Xzz\inz;

% %Edge block
% if nEE > 0
%     inEE=in(nVV+1:nVV+nEE);
%     outEE=pcg(XEE,inEE,1e-8,400);
% else
%     outEE=[];
% end


if nEE > 0
    inEE=in(nVV+1:nVV+nEE);
    tmpE=arg.Lem\(arg.Pem*(inEE));
    tmpE2=arg.Uem\tmpE;
    outEE=arg.Qem*tmpE2;
else
    outEE=[];
end
%Face block

% if nFF > 0
%     inFF=in(nVV+nEE+1:nVV+nEE+nFF);
%     outFF=pcg(XFF,inFF,1e-8,400);
% else
%     outFF=[];
% end


if nFF > 0
    inFF=in(nVV+nEE+1:nVV+nEE+nFF);
    tmpF=arg.Lfm\(arg.Pfm*(inFF));
    tmpF2=arg.Ufm\tmpF;
    outFF=arg.Qfm*tmpF2;
else
    outFF=[];
end


%Interior block 
% if nII > 0
%     inII=in(nVV+nEE+nFF+1:nVV+nEE+nFF+nII);
%     outII=pcg(XII,inII,1e-8,400);
% else
%     outII=[];
% end

if nII > 0
    inII=in(nVV+nEE+nFF+1:nVV+nEE+nFF+nII);
    tmpI=arg.Lim\(arg.Pim*(inII));
    tmpI2=arg.Uim\tmpI;
    outII=arg.Qim*tmpI2;
else
    outII=[];
end


out=[outz;outEg;outFg;outF;outIg;outI;outVV;outEE;outFF;outII];
% out=[outz;outEg;outFg;outF;outIg;outI];
%out=[outVV;outEE;outFF;outII];

%out=[outz;outEg;outFg;outF;outIg;outI];