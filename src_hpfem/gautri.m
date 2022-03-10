function [Quadrature]= gautri(nipd)

% define nip = number of integration points in each direction
%nipd = floor(((2*(order+1))+1+order)/2);

%nipd=ceil((2*(order+1)+max(gorder,order)+1+1)/2);

if nipd<=2
    nipd=2;
end
maxnip = nipd^3;

% obtain integration points......over a tetrahedra
[intxi,inteta,intzeta,intw,nip]=intpoints(nipd,maxnip);
if nip~=nipd^3
    error(message('dimension problem for nipd'));
end

% obtain integration points......over a face
[intfxi,intfet,intfw,nipf]=intpointf(nipd,maxnip);

% Store quadrature data in the Quadrature structure
Quadrature.intxi=intxi;
Quadrature.inteta=inteta;
Quadrature.intzeta=intzeta;
Quadrature.intw=intw;
Quadrature.intfxi=intfxi;
Quadrature.intfet=intfet;
Quadrature.intfw=intfw;
Quadrature.nip=nip;
Quadrature.nipf=nipf;
Quadrature.nipe=nipd;