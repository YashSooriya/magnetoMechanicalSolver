function [IntegratedFields]= PowerAndEnergyComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic,method)

if Options.SplitMech==1 && Options.CustomMRIPost==1
    
    % Get number of output conductivities
    [nCond,~]=size(CondFactorOut);
    % Initialise Output power and disp norm data vectors
    OutPower4K=zeros(nCond*length(freqOut),1);
    OutPower77K=zeros(nCond*length(freqOut),1);
    OutPowerOVC=zeros(nCond*length(freqOut),1);
    DisplacementNorm4K=zeros(nCond*length(freqOut),1);
    DisplacementNorm77K=zeros(nCond*length(freqOut),1);
    DisplacementNormOVC=zeros(nCond*length(freqOut),1);
    
    method = 1;
    
    count=0;
    
    for i=1:nCond
        % Extract conductivities of the specific shields
        condOVC=CondFactorOut(i,1);
        cond77K=CondFactorOut(i,2);
        cond4K=CondFactorOut(i,3);
        
        tic
        
        if method == 1
            [OutPower4K,OutPower77K,OutPowerOVC]=PowerCalculationStaggeredMRIfull(Mesh,Unknown,Basis,Quadrature,ProblemData,UnknownStatic,solStatic,condOVC,cond77K,cond4K,Dynamic,freqOut);
            [DisplacementNorm4K,DisplacementNorm77K,DisplacementNormOVC] = DispNormCalculationMRIfull(Mesh,Unknown,Basis,Quadrature,ProblemData,Dynamic,freqOut);
        else
            for j=1:length(freqOut)
                count=count+1;
                freq=freqOut(j);
                % Extract solution for this conductivity and this frequency for solution matrix
                sol=Dynamic(:,count);
                % Dissipated Power
                [OutPower4K(count),OutPower77K(count),OutPowerOVC(count)]=PowerCalculationStaggeredMRI(Mesh,Unknown,Basis,Quadrature,sol,ProblemData,freq,UnknownStatic,solStatic,condOVC,cond77K,cond4K,Dynamic,freqOut);

                % Calculate displacement norm
                [DisplacementNorm4K(j),DisplacementNorm77K(j),DisplacementNormOVC(j)] = DispNormCalculationMRI(Mesh,Unknown,Basis,Quadrature,sol,ProblemData);

                disp(['------------ Conducting factor choice ',num2str(i),': ',num2str(j),'/',num2str(length(freqOut)), ' frequencies solved. ------------'])

            end
        end
        
        toc
        
    end
    IntegratedFields.OutPower4K=OutPower4K;
    IntegratedFields.OutPower77K=OutPower77K;
    IntegratedFields.OutPowerOVC=OutPowerOVC;
    IntegratedFields.DisplacementNorm4K=DisplacementNorm4K;
    IntegratedFields.DisplacementNorm77K=DisplacementNorm77K;
    IntegratedFields.DisplacementNormOVC=DisplacementNormOVC;
else
    % Get number of output conductivities
    [nCond,~]=size(CondFactorOut);
    % Initialise Output power and disp norm data vectors
    OutPower=zeros(nCond*length(freqOut),ProblemData.NmechBodies);
    DisplacementNorm=zeros(nCond*length(freqOut),ProblemData.NmechBodies);
    count=0;
    for i=1:nCond
        % Extract conductivities of the specific shields
        CondFactorBodies=CondFactorOut(i,:);
        for j=1:length(freqOut)
            count=count+1;
            freq=freqOut(j);
            % Extract solution for this conductivity and this frequency for solution matrix
            sol=Dynamic(:,count);
            % Dissipated Power
                [OutPower(count,:)]=PowerCalculationStaggered(Mesh,Unknown,Basis,Quadrature,sol,ProblemData,freq,UnknownStatic,solStatic,CondFactorBodies);
            % Calculate displacement norm
            [DisplacementNorm(count,:)] = DispNormCalculation(Mesh,Unknown,Basis,Quadrature,sol,ProblemData);
        end
    end
            IntegratedFields.OutPower=OutPower;
        IntegratedFields.DisplacementNorm=DisplacementNorm;
    % Code for only one mechanical problem
end
