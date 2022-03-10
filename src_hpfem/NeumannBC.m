% Function to compute elemental boundary/interface integrals

function [R_A,R_U]= NeumannBC(R_A,R_U,intfw,mue,curle,esize,nm,area,probFlag,probdata,grad_exact,ph,H1bas3D,material,pp)



    
    R_A=R_A+intfw(pp)*(1/(mue))*...
        ( ph(1:esize,1)*(nm(2)*curle(3,1)-nm(3)*curle(2,1)) +...
        ph(1:esize,2)*( nm(3)*curle(1,1)- nm(1)*curle(3,1)) +...
        ph(1:esize,3)*( nm(1)*curle(2,1)-nm(2)*curle(1,1) ))*area;


if probFlag==2
    lambda=probdata.matr.lambda(material);
    G=probdata.matr.G(material);
    div_exact=trace(grad_exact);
    strain_exact=0.5*(grad_exact+grad_exact');
    Idmat=zeros(3,3);
    Idmat(1,1)=1;
    Idmat(2,2)=1;
    Idmat(3,3)=1;
    StressTensor_exact=lambda*div_exact*Idmat+2*G*strain_exact;
    
    
    R_U=R_U-1.256637061435917e-06*intfw(pp)*area*H1bas3D*(StressTensor_exact*nm');
    
end



end



