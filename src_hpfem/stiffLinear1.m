function[K_AA,K_UU]=stiffLinear1(Quadrature,ecphx,ecphy,ecphz,gphx,gphy,gphz,gradH1basisx,gradH1basisy,gradH1basisz,esize,esizeH1,flag,gorder,mue,mycoord,lambda,G,subFlag)

% Extract relevant data from Quadrature structure
intw=Quadrature.intw;
nip=Quadrature.nip;



% intialize stiff matrix
K_AA = zeros(esize,esize);
K_UU = zeros(3*esizeH1,3*esizeH1);
sizeMech=3*esizeH1;

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
        % Commented while solving just the mechanical problem
        
        ph(1:esize,1:3)=(ecphx(i,1:esize)'*asxi(1:3))+...
            (ecphy(i,1:esize)'*aseta(1:3))+...
            (ecphz(i,1:esize)'*aszeta(1:3));
        
        K_AA(1:esize,1:esize)=K_AA(1:esize,1:esize)+...
            ((1/mue)*ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw(i));
        
        if subFlag==2
            expressionforke=zeros(esizeH1,esizeH1,3,3);
        
         gradH1bas(1:esizeH1,1:3)=(gradH1basisx(i,1:esizeH1)'*axi(1:3))+...
            (gradH1basisy(i,1:esizeH1)'*aeta(1:3))+...
            (gradH1basisz(i,1:esizeH1)'*azeta(1:3));
        

        for a=1:esizeH1
            for ii=1:3
                for b=1:esizeH1
                    for jj=1:3
        
        for k=1:3
            for l=1:3
                Csym=Csymfun(ii,jj,k,l,lambda,G);
                expressionforke(a,b,ii,jj)=expressionforke(a,b,ii,jj)+gradH1bas(a,k)*Csym*gradH1bas(b,l);
            end
        end
        localrownumber=(a-1)*3+ii;
        localcolnumber=(b-1)*3+jj;
        
        K_UU(localrownumber,localcolnumber)=K_UU(localrownumber,localcolnumber)+expressionforke(a,b,ii,jj)*det*intw(i);
                    end
                end
            end
        end

       
        
        end
        
        
    end
    
else
    
    for i=1:nip
        %evaluate covairant mapping
        gph(1:gesizet,1)=gphx(i,1:gesizet)';
        gph(1:gesizet,2)=gphy(i,1:gesizet)';
        gph(1:gesizet,3)=gphz(i,1:gesizet)';
        
        % evaluate covairant mapping
            [axi,aeta,azeta,asxi,aseta,aszeta,det]=jacobian_pre(flag,gesizet,gph,mycoord);
            % Commented while solving just the mechanical problem
        
        ph(1:esize,1:3)=(ecphx(i,1:esize)'*asxi(1:3))+...
            (ecphy(i,1:esize)'*aseta(1:3))+...
            (ecphz(i,1:esize)'*aszeta(1:3));
        
        K_AA(1:esize,1:esize)=K_AA(1:esize,1:esize)+...
            ((1/mue)*ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw(i));
        
        if subFlag==2
            expressionforke=zeros(esizeH1,esizeH1,3,3);
               gradH1bas(1:esizeH1,1:3)=(gradH1basisx(i,1:esizeH1)'*axi(1:3))+...
            (gradH1basisy(i,1:esizeH1)'*aeta(1:3))+...
            (gradH1basisz(i,1:esizeH1)'*azeta(1:3));
        
         for a=1:esizeH1
            for ii=1:3
                for b=1:esizeH1
                    for jj=1:3
        
        for k=1:3
            for l=1:3
                Csym=Csymfun(ii,jj,k,l,lambda,G);
                expressionforke(a,b,ii,jj)=expressionforke(a,b,ii,jj)+gradH1bas(a,k)*Csym*gradH1bas(b,l);
            end
        end
        localrownumber=(a-1)*3+ii;
        localcolnumber=(b-1)*3+jj;
        
        K_UU(localrownumber,localcolnumber)=K_UU(localrownumber,localcolnumber)+expressionforke(a,b,ii,jj)*det*intw(i);
                    end
                end
            end
        end
        

        
        end
        
    end
    
    
    
end