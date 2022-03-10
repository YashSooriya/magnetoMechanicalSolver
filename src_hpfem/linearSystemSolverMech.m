function [U]=linearSystemSolverMech(Unknown,Ared,w0,w1,w2,...
    M,K,Resid,dampRatio,NmechBodies)

%=========================================================================
% Initialisations
%=========================================================================
lenEM=Unknown.EM.nunkt;
lenM=Unknown.Mech.nunkt;
nunkC=Unknown.Mech.nunkC;
ndirEM=Unknown.EM.npec;
ndirMech=Unknown.Mech.npec;
%=========================================================================
% System assembly
%=========================================================================

% Construct the Tangent stiffness matrix
TSM = w2*M + w0*K;

wReg=(-1i)*w1;
alpha_M=2*wReg*dampRatio;

clear C
clear CReg


clear Cpre
clear K
clear CpreReg

%=========================================================================
% Re-structure the matrix blocks to solve as fixed point
%=========================================================================
% extract all the blocks: AA, AU, UA, UU from TSM
lenCont=lenEM;
TSM_UA=cell(NmechBodies,1);
for i=1:NmechBodies
TSM_UA{i}=TSM(lenCont+1:lenCont+nunkC(i),1:lenEM);
lenCont=lenCont+nunkC(i);
end
TSM_UU=TSM(lenEM+1:end,lenEM+1:end)+w1*alpha_M*M(lenEM+1:end,lenEM+1:end);
%TSM_UU=TSM(lenEM+1:end,lenEM+1:end);

clear TSM
clear M
tic
%tStartEig=tic;
%EigValues4K=eigs(TSM_UU(1:nunk4K-lenEM,1:nunk4K-lenEM),20,'smallestabs');
%EigValues77K=eigs(TSM_UU(nunk4K-lenEM+1:nunk77K-lenEM,nunk4K-lenEM+1:nunk77K-lenEM),20,'smallestabs');
%EigValuesOVC=eigs(TSM_UU(nunk77K-lenEM+1:nunkOVC-lenEM,nunk77K-lenEM+1:nunkOVC-lenEM),20,'smallestabs');
%TimeEigenAnalysis=toc(tStartEig);
%disp(['Eigenvalue Analysis time = ' num2str(TimeEigenAnalysis)]);
    
%-------------------------------------------------------------------------
% Solve Mechanical problem
%------------------------------------------------------------------------
lenCont=lenEM;
U_UU=[];
for i=1:NmechBodies
    Resid_i=Resid(lenCont+1:lenCont+nunkC(i))+TSM_UA{i}*Ared;
    lenCont=lenCont+nunkC(i);
    U_UUi=TSM_UU(lenCont-lenEM-nunkC(i)+1:lenCont-lenEM,lenCont-lenEM-nunkC(i)+1:lenCont-lenEM)\-Resid_i;             
    U_UU=[U_UU;U_UUi];
end
 
toc
disp('for solving the mechanical problem at this specific frequency.')
TimeSolver=toc;
% Include the Dirichlet Dof into the update vector U, where Dir=0;
U=[Ared(1:lenEM,1);zeros(ndirEM,1);U_UU;zeros(ndirMech,1)];
% U=[U(1:lenEM);zeros(ndirEM,1);U(lenEM+1:end);zeros(ndirMech,1)];


end