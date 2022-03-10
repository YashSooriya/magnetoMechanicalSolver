function[C_AA,C_AU]=dampLinear(C_AA,C_AU,esize,det,ph,intw,curlph,Adc,probstatic,subFlag,coupling,ProblemData,material,H1bas3D)




C_AA(1:esize,1:esize)=C_AA(1:esize,1:esize)+(ph(1:esize,1:3)*ph(1:esize,1:3)'*det*intw);

if coupling==1
    if probstatic==0
        if subFlag==2
            sigma=ProblemData.matr.sigma(material);
            curlADC=curlph'*Adc;
            %==========================================================================
            % Vectorised implementation
            %==========================================================================
            B1=[0 curlADC(3) -curlADC(2); -curlADC(3) 0 curlADC(1); curlADC(2) -curlADC(1) 0];
            C_AU=C_AU-1.256637061435917e-06*sigma*ph*B1*H1bas3D*det*intw;
            %==========================================================================
        end
    end
end

