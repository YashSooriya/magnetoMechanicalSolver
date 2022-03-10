function[X]=initialGuess(mesh,Basis,Quadrature,unknown,ProblemData,omega)

% Extract data from structure
nunkt=unknown.nunkt;
npec=unknown.npec;

X=zeros(nunkt+npec,1);

if ProblemData.probFlag==1
[known]= DirichletEM(mesh,Basis,unknown,Quadrature,1,ProblemData,omega);
elseif ProblemData.probFlag==2
    [known]= DirichletMech(mesh,Basis,unknown,Quadrature,ProblemData,omega);
end


X(nunkt+1:nunkt+npec)=known;


