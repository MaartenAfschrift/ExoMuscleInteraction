%% Find optimal exo control for a more complex model
%---------------------------------------------------

clear all; close all; clc;

%% path information

MainDir     = 'C:\Users\u0088756\Documents\PostDoc\Software\ExoSimulation\BenchMarkTest';
Datapath    = fullfile(MainDir,'Data','Gait92Model');
S.IK_path   = fullfile(Datapath,'KS_Refwalk.mot');
S.ID_path   = fullfile(Datapath,'ID_RefWalk.sto');
S.model_path= fullfile(Datapath,'scaled_model_marker.osim');
load(fullfile(Datapath,'RefWalk.mat'));
S.time      = [RefWalk.Event.RHS_COP(10) RefWalk.Event.RHS_COP(10+1)];     % two heelstrikes
S.OutPath   = fullfile(Datapath,'results');

%% Additional software
% TO Run this example, you have to add the project on solving the muscle
% redundancy problem and metabolic energy consumptions to you rmatlab path

% Solve muscle redundancy (De Groote 2016): https://github.com/antoinefalisse/solvemuscleredundancy_dev
addpath(genpath('C:\Users\u0088756\Documents\GitProjects\solvemuscleredundancy'));
% Computation metabolic energy: https://github.com/MaartenAfschrift/MetabolicEnergy_Simulation
addpath(genpath('C:\Users\u0088756\Documents\GitProjects\MetabolicEnergy\MetabolicEnergy\Models'));

%% Input information
Misc.MuscleNames_Input={};      % Selects all muscles for the Input DOFS when this is left empty.
Misc.DofNames_Input={'ankle_angle_r','knee_angle_r','hip_flexion_r','hip_rotation_r','hip_adduction_r'};    % dof for solving muscl redunancy
S.nDof = length(Misc.DofNames_Input);

Misc.f_cutoff_IK = 6;   Misc.f_cutoff_ID=6;
Misc.f_order_ID=6;      Misc.f_cutoff_lMT=6;
Misc.f_order_lMT=6;     Misc.f_cutoff_dM=6;
Misc.f_order_dM=6;      Misc.f_cutoff_IK=6;
Misc.f_order_IK=6;      Misc.Mesh_Frequency=100;

NMuscles = 43;
MuscleInds_S =32:34;    S.MuscleInds_S=MuscleInds_S;

% Exoskeleton settings
S.N_exo_dof     = 3;                      % only ankle actuation
S.Texo          = [0 0 0.8*72];                % peak exo moment (Nm) (Zhang et al. 2017)
S.MaxTMult      = [ 1 1 0];                % for ankle joint, only plantarflexion assistance, no dorsiflexion
S.MinTMult      = [ -1 -1 -1];
S.MaxPower      = [200 200 200];         % peak exoskeleton power (W)
S.ExoDof        = {'knee_angle_r','hip_flexion_r','ankle_angle_r'};
S.IndExo        = [];

% weights in objective function
auxdata.w0 = 1;      % weight on activations
auxdata.w1 = 100;     % weight on reserve actuators
auxdata.w2 = 0.1;      % weight on vA
auxdata.w3 = 0;%0.01;        % weight on metabolic energy

% musculo-tendon settings
S.ATendon               = ones(1,NMuscles)*35;
S.ATendon(MuscleInds_S) = 10;
S.Topt_res              = 1;

S.CreateAdigatorFiles = 1;      % create new adigator files for the continuous and endpoint function ?
S.RunMuscleAnalysis   = 0;      % Run muscle analysis (again ?)

% settings for metabolic energy (update read this from files using Antoine his function)
S.ST_ratio = [];%[0.54 0.53 0.55 0.5 0.38 0.50 0.56 0.80 0.7];  % the 18m model
S.tension  = [];%[0.62 1    0.75 1.5 0.75 0.50 0.69 0.62 0.75]; % the 18m model                      
S.b_Metab  = 10;

% Loop over tendon stiffness vectors
StiffnessVect = [6 8 10 12 15 17 20 25 30 35 40 70 100 200];    % simulate for different values of tendo stiffness  
nSim=length(StiffnessVect);

% setup NLP
[setup,DatStore] = nExo_Sim_Batch_Setup(S,auxdata,Misc);


%% Run simulations with AFO
disp(' ');  disp('Start Run Simulations');  disp(' ');
for i=1:nSim
    if i>1
        S.CreateAdigatorFiles = 0;  % create new adigator files for the continuous and endpoint function ?
        S.RunMuscleAnalysis = 0;
    end
    setup.auxdata.ATendon(MuscleInds_S)=StiffnessVect(i);
    [res,Texo,TID,MusclePower,lMtilde,MActivation,lTtilde,energy_total,~,...
        ~,~,ExoPosWork,ExoNetWork] = nExo_Sim_Batch_Run(setup,DatStore);
    if i==1
        % pre allocate variables
        nColl=length(res.time);
        Texo_Vect=zeros(nColl,S.N_exo_dof,nSim);
        TID_Vect=zeros(nColl,S.nDof,nSim);
        MusclePower_Vect=zeros(nColl,NMuscles,nSim);
        lMtilde_Vect=zeros(nColl,NMuscles,nSim);
        MActivation_Vect=zeros(nColl,NMuscles,nSim);
        lTtilde_Vect=zeros(nColl,NMuscles,nSim);
        E_Vect=zeros(nColl,NMuscles,nSim);        
        Work_Vect   = zeros(S.N_exo_dof,nSim);
        PosWork_Vect   = zeros(S.N_exo_dof,nSim);
    end    
    Texo_Vect(:,:,i)=Texo;
    TID_Vect(:,:,i)=TID;
    MusclePower_Vect(:,:,i)=MusclePower;
    lMtilde_Vect(:,:,i)=lMtilde;
    MActivation_Vect(:,:,i)=MActivation;
    lTtilde_Vect(:,:,i)=lTtilde;
    E_Vect(:,:,i)=energy_total;   
    PosWork_Vect(:,i)=ExoPosWork;
    Work_Vect(:,i)=ExoNetWork;
end

%% Run simulations without with AFO
S.Texo = [0 0 0];                % peak exo moment (Nm)
[setupNoExo,DatStore] = nExo_Sim_Batch_Setup(S,auxdata,Misc);

for i=1:nSim
    if i>1
        S.CreateAdigatorFiles = 0;  % create new adigator files for the continuous and endpoint function ?
        S.RunMuscleAnalysis = 0;
    end
    
    setupNoExo.auxdata.ATendon(MuscleInds_S)=StiffnessVect(i);
    [res,Texo,TID,MusclePower,lMtilde,MActivation,lTtilde,energy_total,~,...
        ~,~,ExoPosWork,ExoNetWork] = nExo_Sim_Batch_Run(setupNoExo,DatStore);    
    if i==1
        Texo_Vect2=zeros(nColl,S.N_exo_dof,nSim);
        TID_Vect2=zeros(nColl,S.nDof,nSim);
        MusclePower_Vect2=zeros(nColl,NMuscles,nSim);
        lMtilde_Vect2=zeros(nColl,NMuscles,nSim);
        MActivation_Vect2=zeros(nColl,NMuscles,nSim);
        lTtilde_Vect2=zeros(nColl,NMuscles,nSim);
        E_Vect2=zeros(nColl,NMuscles,nSim);
        Work_Vect2   = zeros(S.N_exo_dof,nSim);
        PosWork_Vect2   = zeros(S.N_exo_dof,nSim);
    end
    

    Texo_Vect2(:,:,i)=Texo;
    TID_Vect2(:,:,i)=TID;
    MusclePower_Vect2(:,:,i)=MusclePower;
    lMtilde_Vect2(:,:,i)=lMtilde;
    MActivation_Vect2(:,:,i)=MActivation;
    lTtilde_Vect2(:,:,i)=lTtilde;
    E_Vect2(:,:,i)=energy_total;
    PosWork_Vect2(:,i)=ExoPosWork;
    Work_Vect2(:,i)=ExoNetWork;
end

%% Data export

save(fullfile(S.OutPath,'Energy_TendonSim.mat'),'S','Misc','auxdata','Texo_Vect','TID_Vect','MusclePower_Vect',...
    'lMtilde_Vect','MActivation_Vect','lTtilde_Vect','E_Vect',...
    'Texo_Vect2','TID_Vect2','MusclePower_Vect2','lMtilde_Vect2','MActivation_Vect2',...
    'lTtilde_Vect2','res','DatStore','setup','setupNoExo','E_Vect2','StiffnessVect',...
    'PosWork_Vect2','PosWork_Vect','Work_Vect2','Work_Vect','StiffnessVect');
