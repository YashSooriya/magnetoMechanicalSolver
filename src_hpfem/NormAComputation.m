function [IntegratedFields]= NormAComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic)

if Options.SplitMech==1 && Options.CustomMRIPost==1
    
    % Get number of output conductivities
    [nCond,~]=size(CondFactorOut);
    % Initialise Output power and disp norm data vectors
    OutVectorPotential4K=zeros(nCond*length(freqOut),1);
    OutVectorPotential77K=zeros(nCond*length(freqOut),1);
    OutVectorPotentialOVC=zeros(nCond*length(freqOut),1);
    
    count=0;
    for i=1:nCond
        % Extract conductivities of the specific shields
        condOVC=CondFactorOut(i,1);
        cond77K=CondFactorOut(i,2);
        cond4K=CondFactorOut(i,3);
        for j=1:length(freqOut)
            count=count+1;
            freq=freqOut(j);
            % Extract solution for this conductivitOutPowerOVCy and this frequency for solution matrix
            sol=Dynamic(:,count);
            % Vector Potential
            [OutVectorPotential4K(count),OutVectorPotential77K(count),OutVectorPotentialOVC(count)]=NormACalculationStaggeredMRI(Mesh,Unknown,Basis,Quadrature,sol,ProblemData,freq,UnknownStatic,solStatic,condOVC,cond77K,cond4K);
            disp(['------------ Conducting factor choice ',num2str(i),': ',num2str(j),'/',num2str(length(freqOut)), ' frequencies solved. ------------'])
        end
    end
    IntegratedFields.OutNormA4K=OutVectorPotential4K;
    IntegratedFields.OutNormA77K=OutVectorPotential77K;
    IntegratedFields.OutNormAOVC=OutVectorPotentialOVC;
else
    % Get number of output conductivities
    [nCond,~]=size(CondFactorOut);
    % Initialise Output power and disp norm data vectors
    OutVectorPotential=zeros(nCond*length(freqOut),ProblemData.NmechBodies);
    count=0;
    for i=1:nCond
        % Extract conductivities of the specific shields
        CondFactorBodies=CondFactorOut(i,:);
        for j=1:length(freqOut)
            count=count+1;
            freq=freqOut(j);
            % Extract solution for this conductivity and this frequency for solution matrix
            sol=Dynamic(:,count);
            % Vector Potential
            [OutVectorPotential(count,:)]=NormACalculationStaggered(Mesh,Unknown,Basis,Quadrature,sol,ProblemData,freq,UnknownStatic,solStatic,CondFactorBodies);
            disp(['------------ Conducting factor choice ',num2str(i),': ',num2str(j),'/',num2str(length(freqOut)), ' frequencies solved. ------------'])
        end
    end
    IntegratedFields.OutNormA=OutVectorPotential;
    % Code for only one mechanical problem
end
