function mainParallel(N_s)
 
orderEM = 2;
orderMech = 2;
CondFactorSample = [1 1 1];
CondFactorOut = [1 1 1];
CondFactorChoice = 'default';
chooseSample = "marcos";
chooseOut = "single";
dampRatio = 2e-3;
dampChoice = '2e-3';
nModes = 20;
Ncores = 2;

if chooseSample == "marcos"
    % ns 2324 1296 599 359 180 90 45 23
    % delf1 1 2 5 5 10 20 40 80
    % delf2 3 5 10 25 50 100 200 400
    del_f1 = 10;
    del_f2 = 50;
    s1 = 10:del_f1:1000;
    s2 = 1000:del_f2:5000;
    freqSample = [s1, s2(1:end-1)];
elseif chooseSample == "evenly spaced"
    freqSample = linspace(5, 5000, N_s);
elseif chooseSample == "log spaced"
    freqSample = logspace(log10(5), log10(5000), N_s);
end


if chooseOut == "marcos"
    del_fout = 10;
    N_o = 499;
    freqOut = linspace(15,15+(N_o-1)*del_fout, N_o);
elseif chooseOut == "evenly spaced"
    N_o = 40;
    freqOut = linspace(15, 5000, N_o);
elseif chooseOut == "single"
    freqOut = 1000;
elseif chooseOut == "snapshots"
    freqOut = freqSample;
end

%====================================================================================================================================================
% Main script used to run the coupled solver (parallel version)
%====================================================================================================================================================
% Input parameters:
% orderEM:               Polynomial order of the H(curl) basis functions
% orderMech:             Polynomial order of the H^1 basis functions
% CondFactorSample:      Matrix containing the conductivities to sample (the factors that multiply the reference conductivity)
%                        Each row is a different choice and each column a different conducting subdomain (OVC, 77K, 4K for MRI)
% CondFactorOut:         Matrix containing the conductivities for which to compute the solution (the factors tha multiply the reference conductivity)
% CondFactorChoice:      User specified label for the chosen output conductivities (Used only to create names for data files)
% freqSample:            Vector containing the frequencies in Hz to sample for each conductivity
% freqOut:               Vector containing the frequencies in Hz for which to compute the solution for each conductivity
% dampRatio:             Damping ratio used to damp the mechanica system    
% dampChoice:            User specified label for the chosen damping ratio (Used only to create names for data files)
% nModes:                Number of POD modes
% Ncores:                Number of workers to use for parallel computation
%======================================================================================================================================================
clc
close all
format long e
addpath(genpath('./'))

%-------------------------------------------------------------------------
% Problem Switches
%-------------------------------------------------------------------------
% Switches to turn certain feature on/off
% 1 - on, 0 - off

% Pre processing flags
ReadMesh         = 1;          % Read mesh (1) or load existing mesh data (0)

% Solver Flags
POD              = 0;          % Use POD ROM (1) or Full Order (0)
PODP             = 0;          % Use PODP (1) or PODI(0). Note PODI has been implemented for 1 parameter only.
Offline          = 1;          % Compute off-line POD stage (1) or load existing data from file (0)
SourceMapping    = 0;          % Map current to solenoidal space (1) or not (0) (Required for transversal coils)
Assemble         = 1;          % Assemble system matrices (1) or load assembled matrices from file
couple           = 0;          % Solve coupled (1) or decoupled (0) problem
SplitMech        = 0;          % Split the mechanical problem in three separated problems for the OVC, 77K and 4K (customised for these problems)
StaticMechanics  = 1;          % Run (1) or not (0) the static solver for the mechanical field (Only on if magnetic material)
Non0Dir          = 1;          % Consider non-zero Dirichlet values (1) or not (0)
freqSweep        = 1;          % Assemble Dirichlet DOF using frequency sweep method (requires a lot of memory for large problems)
                               % Only necessary for non zero Dirichlet BC
% Post processing flags
offlineOnly      = 0;          % Only calculate the offline stage, skip the online (only for POD ROM)
fieldCalc        = 0;          % Calculate integrated field quantities (Output Power, Kinetic Energy)
CustomMRIPost    = 0;          % (1)Return dissipated power and kinetic energy directly for 4K, 77K and OVC shields (Customised for MRI problems)
                               % (0) Return the dissipaed power and kinetic energy for the different unnamed mechanical subdomains

normA            = 0;          % Return the integrated Magnetic vector potential for 4K, 77K and OVC shields (1) or not (0)
normCurlA        = 0;          % Return the integrated curl of the Magnetic vector potential for 4K, 77K and OVC shields (1) or not (0)
linePlotOn       = 0;          % Line plot of the fields
paraview         = 1;          % Paraview .vtk file writer
errorsOn         = 0;          % Compute error with respect to analytical solution
ErrorPOD         = 0;          % Compute error of POD with respect to full order (e_2)
FixedPoint       = 0;

% svd save file name
svdSave          = 0;          % Use custom save name for SVD, will use default if 0
svdSaveName      = 'SVDResult_Ns40_m20_105hz.mat';% File name for SVD mat file

%-------------------------------------------------------------------------
% Problem Definition
%-------------------------------------------------------------------------

% Define string with problem file name
%problem= 'MHIGradXSplitNewCoil';
problem='Cylinder';

%=========================================================================
% Extract the problem data from the problem file
%=========================================================================

% Determine problem data based on the problem files
ProblemData = eval(['problem',num2str(problem),'(orderEM)']);

% Store the extra problem specific data in the ProblemData structure
ProblemData.order = orderEM;
ProblemData.orderH1 = orderMech;
MaxOrder=max(orderEM,ProblemData.orderH1);

%-------------------------------------------------------------------------
% Read mesh
%-------------------------------------------------------------------------

% Job data
job = ProblemData.jb.job;
fileName=strcat('Mesh_',job);

if ReadMesh==1
if ProblemData.jb.meshtype==3
    meshdata = textread([job '.vol'], '%s','delimiter', '\n');
    disp('reading netgen mesh....')
    [Mesh]=meshinfo(meshdata,ProblemData);
    disp('...done')
elseif ProblemData.jb.meshtype==4
    disp('reading ANSYS mesh....')
    [Mesh]= ansysmeshinfo(ProblemData);
    disp('...done')
elseif ProblemData.jb.meshtype==5
    name=ProblemData.jb.name;
    disp('Reading Opera mesh...')
    [Mesh] = OperaMeshReader(name,ProblemData);
    disp('...done')
end
    save(fileName,'Mesh');
else
    load (fileName,'Mesh');
end



% Flag those elements in conducting regions to include gradients of basis
[Mesh]=ConductorBodiesFlag(ProblemData,Mesh);

% Define size of elemental arrays for H(curl) and H^1 spaces
[ProblemData]=GetSizeArrays(ProblemData,Mesh);

%-------------------------------------------------------------------------
% Extra mesh data
%-------------------------------------------------------------------------
% Set reference element types (Type I or Type II)
[Mesh]=gettype(Mesh);
disp('Completed setting reference types')

% Compute the numbering of the edges
[Mesh]=edgeno(Mesh);
disp('Completed edge numbering')

% Compute numbering of faces and asign boundary conditions 
[Mesh]=AssignBC(Mesh,ProblemData);
disp('Completed face numbering')

% Extract local mesh for the mechanical problem
[Mesh.mech]=localMesh(Mesh,ProblemData,2);

%========================================================================
% Compute properties of the numerical integration in reference element
%========================================================================
% Compute the quadrature points, basis functions and blending functions on
% the reference Iso-Parametric element

%--------------------------------------------------------------------------
% Blending function
%--------------------------------------------------------------------------
[Mesh]=getcoeff(Mesh,ProblemData);

display(['Program using orders q = ',num2str(ProblemData.order),' p = ',num2str(ProblemData.orderH1)]);
display(['Size of elemental arrays = ',num2str(ProblemData.esizet)]);

%--------------------------------------------------------------------------
% Gaussian Quadrature
%--------------------------------------------------------------------------
% Define number of integration points
nip=ceil((2*MaxOrder+1)/2);

% Compute the integration weigths and locations on reference element
[Quadrature]=gautri(nip);
%--------------------------------------------------------------------------
% Basis Functions
%--------------------------------------------------------------------------
% Evaluate basis functions at integration points
[Basis]= evalBasis(Quadrature,ProblemData);

%==========================================================================
% Compute the system unknowns of the individual physics types
%==========================================================================
[Unknown]=elementUnknownNumbering(Mesh,ProblemData,SplitMech);

%==========================================================================
% Store program options in the options data structure
%==========================================================================
Options.noSnapshots=size(freqSample,2);
Options.orderEM=orderEM;
Options.orderMech=orderMech;
Options.CondFactorSample=CondFactorSample;
Options.CondFactorOut=CondFactorOut;
Options.CondFactorChoice=CondFactorChoice;
Options.chooseSample=chooseSample;
Options.chooseOut=chooseOut;
Options.dampRatio=dampRatio;
Options.dampChoice=dampChoice;
Options.nModes=nModes;
Options.Couple=couple;
Options.StaticMechanics=StaticMechanics;
Options.fieldCalc=fieldCalc;
Options.linePlotOn=linePlotOn;
Options.Paraview=paraview;
Options.errorsOn=errorsOn;
Options.Assemble=Assemble;
Options.POD=POD;
Options.Offline=Offline;
Options.offlineOnly=offlineOnly;
Options.SourceMapping=SourceMapping;
Options.SplitMech=SplitMech;
Options.CustomMRIPost=CustomMRIPost;
Options.freqSweep=freqSweep;
Options.Non0Dir=Non0Dir;
Options.normA=normA;
Options.svdSave=svdSave;
Options.svdSaveName=svdSaveName;
%=======================================================================================================================================
% Coupled Problem solver
%=======================================================================================================================================
% Start the timer
tic

% Static Solver
[Static,UnknownStatic]=staticSolver(Mesh,Basis,Quadrature,ProblemData,Options);
% Source mapping to compute JAC in terms of a vector potential
if SourceMapping==1
[StaticCurrent,UnknownCurrent]=staticSolver2(Mesh,Basis,Quadrature,ProblemData,Options,Static);
else
    StaticCurrent=Static;
    UnknownCurrent=UnknownStatic;
end

% Dynamic Solver (Frequency domain)
if POD==1
    if PODP==1
        [Dynamic]=frequencySolverPODPParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,dampChoice,dampRatio,CondFactorOut,CondFactorChoice,CondFactorSample,freqSample,nModes,Ncores);
    else
        [Dynamic]=frequencySolverPODIParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,dampChoice,dampRatio,CondFactorOut,CondFactorChoice,CondFactorSample,freqSample,nModes,Ncores);
    end
else
    [Dynamic]=frequencySolverFullParallel(Static,StaticCurrent,UnknownCurrent,UnknownStatic,Mesh,Basis,Quadrature,Unknown,ProblemData,Options,freqOut,dampChoice,dampRatio,CondFactorOut,CondFactorChoice,Ncores);
end


%==========================================================================
% Post processing
%==========================================================================
clear StaticCurrent UnknownCurrent
% Extract static solution
solStatic=Static.sol(:,1);

if  linePlotOn==1
    % Line Plots
    myplot3(Mesh,Unknown,ProblemData,Dynamic(:,1),problem,0,freqRange)
end

%=========================================================================
% ParaView (3D plots)
%=========================================================================
if paraview==1
    % 3D plots for the dynamic solution (Paraview)
    if FixedPoint==1
        pointvalue_multi(Mesh,Unknown,ProblemData,Dynamic(:,1),0,freqRange);
    else
        pointvalue_multi_Staggered(Mesh,Unknown,ProblemData,Dynamic(:,1),0,freqOut,UnknownStatic,solStatic);
      %  pointvalue_multi_Staggered_Air(Mesh,Unknown,ProblemData,Dynamic(:,1),0,freqOut);
    end
    % 3D plots for the static solution (Paraview)
    pointvalue_multi(Mesh,UnknownStatic,ProblemData,solStatic,1,freqOut);
end
%==========================================================================




%=========================================================================
% Compute error of solution (Full order) if analytical solution is available
%=========================================================================
if errorsOn==1
     hcurlerr=zeros(length(freqOut),1);
     [L2_error_mech,Energy_error_mech,SNS_error_mech] =ErrorStatic(Mesh,UnknownStatic,Basis,Quadrature,solStatic,ProblemData,1,50);
    for j=1:length(freqOut)
        freq=freqOut(j);
        sol=Dynamic(:,j);
        [hcurlerr(j),DisplacementNorm(j)] =ErrorDynamic(Mesh,Unknown,Basis,Quadrature,sol,ProblemData,0,freq);
    end
end
%=========================================================================

%=========================================================================
% Compute dissipated power and kinetic Energy
%=========================================================================
if fieldCalc==1
    disp('Calcuating Power and Energy ...')
    [IntegratedFields]= PowerAndEnergyComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic);
    currDate = strrep(datestr(datetime), ':','_');
    folder = ['data/powerEnergy/',currDate,'/'];
    mkdir(folder)
    writestruct(Options,[folder,'Options.xml'])
    saveFile=[folder,'FrequencySweepMHIGradXPowerEnergy'];
    save(saveFile,'IntegratedFields');
    disp(['Saved to ', saveFile])
end

%=========================================================================
% Compute the Norm of A field
%=========================================================================

if normA==1
    disp('Calcuating Norm A ...')
    [IntegratedNormA]= NormAComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic);
    currDate = strrep(datestr(datetime), ':','_');
    folder = ['data/normA/',currDate,'/'];
    mkdir(folder)
    writestruct(Options,[folder,'Options.xml'])
    saveFile=[folder,'FrequencySweepMHIGradXNormA'];
    save(saveFile,'IntegratedNormA');
    disp(['Saved to ', saveFile])
end
%=========================================================================

%=========================================================================
% Compute the Norm of the curlA field
%=========================================================================

if normCurlA==1
    disp('Calcuating Norm curlA ...')
    [IntegratedNormCurlA]= CurlNormAComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic);
    currDate = strrep(datestr(datetime), ':','_');
    folder = ['data/normCurlA/',currDate,'/'];
    mkdir(folder)
    writestruct(Options,[folder,'Options.xml'])
    saveFile=[folder,'FrequencySweepMHIGradXNormCurlA'];
    save(saveFile,'IntegratedNormCurlA');
    disp(['Saved to ', saveFile])
end
%=========================================================================

%=========================================================================
% Compute magnetic vector potential
%=========================================================================
if normA==1
    [IntegratedNormA]= NormAComputation(Options,CondFactorOut,freqOut,Mesh,Unknown,Basis,Quadrature,Dynamic,ProblemData,UnknownStatic,solStatic);
    fileName=sprintf('FrequencySweepParallelFullTestNormA%dto%dHz_damp%s_q%d_p%d_CondFactor%sOVCNS',freqOut(1),freqOut(length(freqOut)),dampChoice,orderEM,orderMech,CondFactorChoice);
    save(fileName,'IntegratedNormA');
end

%=========================================================================
% Compute error of POD solution with respect to full order (e_2)
%=========================================================================
if ErrorPOD==1
DynamicPOD=Dynamic;
fileName=sprintf('FullOrderSolution');
load(fileName,'Dynamic');
[nCond,~]=size(CondFactorOut);
lenEM=Unknown.EM.nunkt;
lenM=Unknown.Mech.nunkt;
ndirEM=Unknown.EM.npec;
ndirMech=Unknown.Mech.npec;
nTotal=ndirEM+ndirMech+lenEM+lenM;
ErrorPOD=zeros(nTotal,nCond*length(freqOut));
ErrorPODEM=zeros(lenEM,nCond*length(freqOut));
ErrorPODMech=zeros(lenM,nCond*length(freqOut));
% Set counter to zero
count=0;
for i=1:nCond
    for j=1:length(freqOut)
        count=count+1;
        ErrorPOD(:,count)=norm(DynamicPOD(:,count)-Dynamic(:,count))/norm(DynamicPOD(:,count));
        ErrorPODEM(:,count)=norm(DynamicPOD(1:lenEM,count)-Dynamic(1:lenEM,count))/norm(DynamicPOD(1:lenEM,count));
        ErrorPODMech(:,count)=norm(DynamicPOD(lenEM+ndirEM+1:lenEM+ndirEM+lenM,count)-Dynamic(lenEM+ndirEM+1:lenEM+ndirEM+lenM,count))/norm(DynamicPOD(lenEM+ndirEM+1:lenEM+ndirEM+lenM,count));
    end
end
fileName=sprintf('ErrorPODCase1');
save(fileName,'ErrorPOD','ErrorPODEM','ErrorPODMech');
end
        


  toc
