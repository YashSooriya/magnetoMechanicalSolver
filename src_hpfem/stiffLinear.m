function[K_AA,K_AU,K_UA,K_UU]=stiffLinear(K_AA,K_AU,K_UA,K_UU,esize,mue,D,subFlag,Adc,det,ph,B,intw,constantB,A2,Grad3D,S2,coupling,FlagVolume)



        K_AA(1:esize,1:esize)=K_AA(1:esize,1:esize)+...
            ((1/(mue))*ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw);
        
        if subFlag==2
 
        K_UU(1:end,1:end)= K_UU(1:end,1:end)+1.256637061435917e-06*B'*D*B*det*intw;

        end
        if coupling==1
            if FlagVolume==1
                curlADC=ph'*Adc;
                
                AuxB=ph*curlADC;
                
                constantB(1:3:end,1)=AuxB;
                constantB(2:3:end,2)=AuxB;
                constantB(3:3:end,3)=AuxB;
                
                
                A2(1:3:end,1)=ph(:,1);
                A2(1:3:end,2)=ph(:,2);
                A2(1:3:end,3)=ph(:,3);
                A2(2:3:end,4)=ph(:,1);
                A2(2:3:end,5)=ph(:,2);
                A2(2:3:end,6)=ph(:,3);
                A2(3:3:end,7)=ph(:,1);
                A2(3:3:end,8)=ph(:,2);
                A2(3:3:end,9)=ph(:,3);
                
                
                
                B2=[curlADC(1) 0 0; 0 curlADC(1) 0; 0 0 curlADC(1); curlADC(2) 0 0; 0 curlADC(2) 0; 0 0 curlADC(2); curlADC(3) 0 0; 0 curlADC(3) 0; 0 0 curlADC(3)];
                B3=[curlADC(1) curlADC(2) curlADC(3); 0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3);0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3)];
                
                S=(1/mue)*(A2*(B2+B3)-constantB);
                
                S2(1,:)=S(1:3:end,1);
                S2(2,:)=S(1:3:end,2);
                S2(3,:)=S(1:3:end,3);
                S2(4,:)=S(2:3:end,1);
                S2(5,:)=S(2:3:end,2);
                S2(6,:)=S(2:3:end,3);
                S2(7,:)=S(3:3:end,1);
                S2(8,:)=S(3:3:end,2);
                S2(9,:)=S(3:3:end,3);
                
                
                
               
                K_UA=K_UA+Grad3D*S2*det*intw;
              
                
             end
        end
        
        

    
    
