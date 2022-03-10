function[M_UU]=massLinear(M_UU,subFlag,probdata,material,det,H1bas3D,intw)
 
        

        
        if subFlag==2
            
        rho=probdata.matr.rho(material);

        M_UU(1:end,1:end)= M_UU(1:end,1:end)+1.256637061435917e-06*(H1bas3D'*H1bas3D)*rho*det*intw;
        
        end
    
    
