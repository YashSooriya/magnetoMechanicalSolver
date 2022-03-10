% Function to compute elemental curl curl stiffness matrix

function stiff= stiffmatrix(esize,intw,nip,...
    ecphx,ecphy,ecphz,flag,gorder,mue,gphx,gphy,gphz,mycoord)


% intialize stiff matrix
stiff = zeros(esize,esize);
ph=zeros(esize,3);
gesizet=(gorder+1+1)*(gorder+1+2)*(gorder+1+3)/6;
gph=zeros(gesizet,3);
 if flag==0
     gesizet=4;
 end


% evaluate covairant mapping (for linear geometry mapping is constant)
if flag==0
    gph(1:gesizet,1)=gphx(1,1:gesizet)';
    gph(1:gesizet,2)=gphy(1,1:gesizet)';
    gph(1:gesizet,3)=gphz(1,1:gesizet)';

    [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);

%     
    for i=1:nip
        
        ph(1:esize,1:3)=(ecphx(i,1:esize)'*asxi(1:3))+...
            (ecphy(i,1:esize)'*aseta(1:3))+...
            (ecphz(i,1:esize)'*aszeta(1:3));
        
        stiff(1:esize,1:esize)=stiff(1:esize,1:esize)+...
            ((1/mue)*ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw(i));
        
        
    end
    
    
else
    
    for i=1:nip
        %evaluate covairant mapping
        gph(1:gesizet,1)=gphx(i,1:gesizet)';
        gph(1:gesizet,2)=gphy(i,1:gesizet)';
        gph(1:gesizet,3)=gphz(i,1:gesizet)';
        
        % evaluate covairant mapping
            [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
        
        ph(1:esize,1:3)=(ecphx(i,1:esize)'*asxi(1:3))+...
            (ecphy(i,1:esize)'*aseta(1:3))+...
            (ecphz(i,1:esize)'*aszeta(1:3));
        
        stiff(1:esize,1:esize)=stiff(1:esize,1:esize)+...
            ((1/mue)*ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw(i));
        
        
    end
    
    
    
end