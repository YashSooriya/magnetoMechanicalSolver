function [U]=linearSystemSolverMechNS(Unknown,Ared,w0,w1,w2,...
    M,K,Resid,dampRatio)

%=========================================================================
% Initialisations
%=========================================================================
lenEM=Unknown.EM.nunkt;
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
TSM_UA=TSM(lenEM+1:end,1:lenEM);
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
Resid_U=Resid(lenEM+1:end)+TSM_UA*Ared;
U_UU=TSM_UU\-Resid_U; 
 
toc
TimeSolver=toc
% Include the Dirichlet Dof into the update vector U, where Dir=0;
U=[Ared(1:lenEM,1);zeros(ndirEM,1);U_UU;zeros(ndirMech,1)];
% U=[U(1:lenEM);zeros(ndirEM,1);U(lenEM+1:end);zeros(ndirMech,1)];


end