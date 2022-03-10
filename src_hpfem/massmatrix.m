% Function to compute elemental mass matrix

function mass= massmatrix(esize,intw,nip,...
    ephx,ephy,ephz,flag,gorder,gphx,gphy,gphz,mycoord)


% intialize mass matrix
mass = zeros(esize,esize);
ph=zeros(esize,3);
gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
gph=zeros(gesizet,3);
 if flag==0
     gesizet=4;
 end



if flag==0
% evaluate covairant mapping (for linear geometry mapping is constant)
    gph(1:gesizet,1)=gphx(1,1:gesizet)';
    gph(1:gesizet,2)=gphy(1,1:gesizet)';
    gph(1:gesizet,3)=gphz(1,1:gesizet)';

    [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
    
    for i=1:nip
        % use stored functions
        ph(1:esize,1:3)=(ephx(i,1:esize)'*axi(1:3))+...
            (ephy(i,1:esize)'*aeta(1:3))+...
            (ephz(i,1:esize)'*azeta(1:3));
        
        mass(1:esize,1:esize)=mass(1:esize,1:esize)+...
            (ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw(i));
    end
else
    % recompute the jacobian
    for i=1:nip
        gph(1:gesizet,1)=gphx(i,1:gesizet)';
        gph(1:gesizet,2)=gphy(i,1:gesizet)';
        gph(1:gesizet,3)=gphz(i,1:gesizet)';
        
        % evaluate covairant mapping
    [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
        
        % use stored functions
        ph(1:esize,1:3)=(ephx(i,1:esize)'*axi(1:3))+...
            (ephy(i,1:esize)'*aeta(1:3))+...
            (ephz(i,1:esize)'*azeta(1:3));
        
        mass(1:esize,1:esize)=mass(1:esize,1:esize)+...
            (ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw(i));
    end
end