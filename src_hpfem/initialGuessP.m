function[X]=initialGuessP(mesh,Basis,Quadrature,unknown,ProblemData,omega,probFlag)

% Extract data from structure
nunkt=unknown.nunkt;
npec=unknown.npec;

X=zeros(nunkt+npec,1);

% if probFlag==1
% [known]= DirichletEMPOD(mesh,Basis,unknown,Quadrature,1,ProblemData,omega,probFlag);
% elseif probFlag==2
%     [known]= DirichletMech(mesh,Basis,unknown,Quadrature,ProblemData,omega,probFlag);
% end

if probFlag==1
[known]= DirichletEM(mesh,Basis,unknown,Quadrature,1,ProblemData,omega);
elseif probFlag==2
    [known]= DirichletMech(mesh,Basis,unknown,Quadrature,ProblemData,omega);
end

X(nunkt+1:nunkt+npec)=known;


