% Function to compute the source term inside the conductor

function [erhs,mdiopole]= jmatrix(esize,intw,nip,...
    ephx,ephy,ephz,flag,gorder,ljsrc,erhs,mdiopole,probdata,...
    gphx,gphy,gphz,kappa,mycoord,phh1,domain)

muz=probdata.matr.muz;
nrhs=probdata.jb.nrhs;
gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
gph=zeros(gesizet,3);
ph1=zeros(gesizet,1);
 if flag==0
     gesizet=4;
 end 
ph=zeros(esize,3);


if flag==0
%   evaluate covairant mapping (for linear geometry mapping is constant)

    gph(1:gesizet,1)=gphx(1,1:gesizet)';
    gph(1:gesizet,2)=gphy(1,1:gesizet)';
    gph(1:gesizet,3)=gphz(1,1:gesizet)';
    
    [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);

    for i=1:nip
        ph1(1:gesizet,1)=phh1(i,1:gesizet)';

        % computx x,y,z
        [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
        
        % use stored functions
        ph(1:esize,1:3)=(ephx(i,1:esize)'*axi(1:3))+...
            (ephy(i,1:esize)'*aeta(1:3))+...
            (ephz(i,1:esize)'*azeta(1:3));
        
         % Use problem file to define src function
        fun=probdata.es.srcfun;
        arg=probdata.es.srcfunarg;
        arg.kappa=kappa;
        src=fun(x,y,z,ljsrc,arg,domain);
        
        erhs(1:esize,1:nrhs)=erhs(1:esize,1:nrhs)+...
            (ph(1:esize,1:3)*src(1:3,1:nrhs)*det*intw(i));

    
    end
else
    % recompute the jacobian
    for i=1:nip
        gph(1:gesizet,1)=gphx(i,1:gesizet)';
        gph(1:gesizet,2)=gphy(i,1:gesizet)';
        gph(1:gesizet,3)=gphz(i,1:gesizet)';
        ph1(1:gesizet,1)=phh1(i,1:gesizet)';
        
        % evaluate covairant mapping
        [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord) ;


        % computx x,y,z
        [x,y,z]= getxyzcu_pre(ph1,mycoord,gesizet);
        
        % use stored functions
        ph(1:esize,1:3)=(ephx(i,1:esize)'*axi(1:3))+...
            (ephy(i,1:esize)'*aeta(1:3))+...
            (ephz(i,1:esize)'*azeta(1:3));
        
        % Use problem file to define src function
        fun=probdata.es.srcfun;
        arg=probdata.es.srcfunarg;
        arg.kappa=kappa;
        src=fun(x,y,z,ljsrc,arg,domain);

        
        erhs(1:esize,1:nrhs)=erhs(1:esize,1:nrhs)+...
            (ph(1:esize,1:3)*src(1:3,1:nrhs)*det*intw(i));

    
    end   
end

