% Function to compute elemental boundary/interface integrals

function [K_UA,K_UA2,K_UA3,R_U,StressIntegral]= BoundTermsIntegration(K_UA,K_UA2,K_UA3,R_U,esize,esizeH1,intfw,mue,...
    contador,area,Nhf,nm,pp,constantB,A2,Basis3Df,ConstantD,sigmae,probstatic,StressIntegral,ecph,A,curlADC)




% Compute curl(A) for stress tensor computation

curlADC=ecph'*A;


AuxB=ecph*curlADC;

% Define constant for vectorized definition
constantB(1:3:end,1)=AuxB;
constantB(2:3:end,2)=AuxB;
constantB(3:3:end,3)=AuxB;

% Define the components that we need
A2(1:3:end,1)=ecph(:,1);
A2(1:3:end,2)=ecph(:,2);
A2(1:3:end,3)=ecph(:,3);
A2(2:3:end,4)=ecph(:,1);
A2(2:3:end,5)=ecph(:,2);
A2(2:3:end,6)=ecph(:,3);
A2(3:3:end,7)=ecph(:,1);
A2(3:3:end,8)=ecph(:,2);
A2(3:3:end,9)=ecph(:,3);


% Define matrices for vectorised calculus of coupling
% blocks
 B2=[curlADC(1) 0 0; 0 curlADC(1) 0; 0 0 curlADC(1); curlADC(2) 0 0; 0 curlADC(2) 0; 0 0 curlADC(2); curlADC(3) 0 0; 0 curlADC(3) 0; 0 0 curlADC(3)];
 B3=[curlADC(1) curlADC(2) curlADC(3); 0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3);0 0 0; 0 0 0; 0 0 0; curlADC(1) curlADC(2) curlADC(3)];

% Define linearised Maxwell stress tensor
 S=(1/mue)*(A2*(B2+B3)-constantB);
% Multiply by normal vector
constantD=S*nm';
% Rearrange for vectorised calculus
constantD2=reshape(constantD,3,esize);
constantD3=repmat(constantD2,esizeH1,1);
v = ceil([1:(esizeH1*3)]./3);
helpMatrix=diag(Nhf(v),0);

% Terms for the residual boundary term
Constant=0.5*(curlADC'*curlADC);

% Define full Maxwell stress tensor (necessary components)
sigmae(1)=curlADC(1)*curlADC(1)-Constant;
sigmae(2)=curlADC(2)*curlADC(2)-Constant;
sigmae(3)=curlADC(3)*curlADC(3)-Constant;
sigmae(4)=curlADC(1)*curlADC(2);
sigmae(5)=curlADC(1)*curlADC(3);
sigmae(6)=curlADC(2)*curlADC(3);
sigmae=(1/(mue))*sigmae;

% Multiply by normal vector
ConstantD(1)=sigmae(1)*nm(1)+sigmae(4)*nm(2)+sigmae(5)*nm(3);
ConstantD(2)=sigmae(4)*nm(1)+sigmae(2)*nm(2)+sigmae(6)*nm(3);
ConstantD(3)=sigmae(5)*nm(1)+sigmae(6)*nm(2)+sigmae(3)*nm(3);



% Add coupling contributions to stiffness matrix
% if contador==2
%     
%     K_UA=K_UA-helpMatrix*constantD3*area*intfw(pp);
%     
% elseif contador==1
%     
%     K_UA2=K_UA2-helpMatrix*constantD3*area*intfw(pp);
%     
% elseif contador==3
%     
%     K_UA3=K_UA3-helpMatrix*constantD3*area*intfw(pp);
%     
% end



% Add coupling contributions to residual
if probstatic==1
    %R_U=R_U-Basis3Df*ConstantD*intfw(pp)*area;
    StressIntegral=StressIntegral+ConstantD*intfw(pp)*area;
end


end



