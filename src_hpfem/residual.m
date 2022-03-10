function[R_A,R_U]=residual(R_A,R_U,mue,Adc,A0,probstatic,det,ph2,ph,H1bas3D,intw,current,bodyForce,Grad3D,sigmae,coupling,FlagVolume,SourceMapping)

% Current Source term (EM)
if probstatic==1
    R_A = R_A-det*intw*ph2*current;
elseif probstatic==2
    R_A = R_A-det*intw*ph*current;
else
    if SourceMapping==1
        A0Vector=ph2'*A0;
        R_A = R_A-det*intw*ph*A0Vector;
    else
        R_A = R_A-det*intw*ph2*current;
    end
end


% Body force term (mechanics)
R_U=R_U+det*intw*H1bas3D'*bodyForce*1.256637061435917e-06;


if probstatic==1
    
    if coupling==1
        if FlagVolume==1
            curlADC=ph'*Adc;
            
            Constant=0.5*(curlADC'*curlADC);
            
            sigmae(1)=curlADC(1)*curlADC(1)-Constant;
            sigmae(5)=curlADC(2)*curlADC(2)-Constant;
            sigmae(9)=curlADC(3)*curlADC(3)-Constant;
            sigmae(2)=curlADC(1)*curlADC(2);
            sigmae(3)=curlADC(1)*curlADC(3);
            sigmae(4)=sigmae(2);
            sigmae(6)=curlADC(2)*curlADC(3);
            sigmae(7)=sigmae(3);
            sigmae(8)=sigmae(6);
            sigmae=(1/(mue))*sigmae;
            
            R_U=R_U+Grad3D*sigmae*intw*det;
        end
        
        
    end
    
    
end







