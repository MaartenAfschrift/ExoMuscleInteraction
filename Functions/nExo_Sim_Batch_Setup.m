function [setup,DatStore,auxdata] = nExo_Sim_Batch_Setup(S,auxdata,Misc)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% scale optimal control exo
auxdata.Texo    = S.Texo;       % max acutation = 50 Nm
auxdata.nExo    = S.N_exo_dof;

auxdata.ST_ratio = S.ST_ratio ;  % the 18m model
auxdata.tension  = S.tension ; % the 18m model         

auxdata.b_Metab = S.b_Metab;
auxdata.Topt_res= S.Topt_res;

auxdata.ID_rel = 0;
if isfield(S,'ID_rel')
   auxdata.ID_rel=S.ID_rel;
else
    auxdata.ID_rel=zeros(size(S.Texo));
end

% path information
% PathMuscleRedundancy = 'C:\Users\u0088756\Documents\GitProjects\solvemuscleredundancy';
% addpath(genpath(PathMuscleRedundancy));
% MetabolicEnergyPath = 'C:\Users\u0088756\Documents\GitProjects\MetabolicEnergy\MetabolicEnergy\Models';
% addpath(genpath(MetabolicEnergyPath));
% addpath(genpath(pwd));

MainDir = 'C:\Users\u0088756\Documents\GitProjects\SimulateExo\OptimalControl';
OutData = fullfile(MainDir,'Data','OsimExample');
if ~isdir(OutData); mkdir(OutData); end;

% Needed Input Arguments
IK_path=S.IK_path;
ID_path=S.ID_path;
model_path=S.model_path;
time=S.time;    % Part of the right stance phase
OutPath=S.OutPath;


%% Get motion input information

% muscle analysis
Misc.time=time;
MuscleAnalysisPath=fullfile(OutPath,'MuscleAnalysis'); if ~exist(MuscleAnalysisPath,'dir'); mkdir(MuscleAnalysisPath); end
disp('MuscleAnalysis Running .....');
if S.RunMuscleAnalysis
    OpenSim_Muscle_Analysis(IK_path,model_path,MuscleAnalysisPath,[time(1) time(end)])
end
disp('MuscleAnalysis Finished');
Misc.MuscleAnalysisPath=MuscleAnalysisPath;

% ----------------------------------------------------------------------- %
% Extract muscle information -------------------------------------------- %

% Get number of degrees of freedom (dofs), muscle-tendon lengths and moment
% arms for the selected muscles.
[~,Misc.trialName,~]=fileparts(IK_path);
if ~isfield(Misc,'MuscleNames_Input') || isempty(Misc.MuscleNames_Input)
    Misc=getMuscles4DOFS(Misc);
end
[DatStore] = getMuscleInfo(IK_path,ID_path,Misc);

Misc.ATendon = S.ATendon;

% get specific tension and fiber type composition
%------------------------------------------------
if isempty(auxdata.ST_ratio)
    ST_ratio = getSlowTwitchRatios(DatStore.MuscleNames); 
    auxdata.ST_ratio = ST_ratio';
end
if isempty(auxdata.tension)
    tension = getSpecificTensions(DatStore.MuscleNames);
    auxdata.tension = tension';
end
if isfield(S,'ExoDof')
    auxdata.IndExo=ones(1,length(S.ExoDof));
    for i=1:length(S.ExoDof)
        auxdata.IndExo(i) =  find(strcmp(S.ExoDof{i},DatStore.DOFNames));
    end
end
      

% Extract the muscle-tendon properties
[DatStore.params,DatStore.lOpt,DatStore.L_TendonSlack,DatStore.Fiso,DatStore.PennationAngle]=ReadMuscleParameters(model_path,DatStore.MuscleNames);
% Static optimization using IPOPT solver
DatStore = SolveStaticOptimization_IPOPT(DatStore);

%% Setup optimal control problem

% Input arguments
auxdata.NMuscles = DatStore.nMuscles;   % number of muscles
auxdata.Ndof = DatStore.nDOF;           % number of dofs
auxdata.ID = DatStore.T_exp;            % inverse dynamics
auxdata.params = DatStore.params;       % Muscle-tendon parameters
auxdata.IK = DatStore.q_exp;

% ADiGator works with 2D: convert 3D arrays to 2D structure (moment arms)
for i = 1:auxdata.Ndof
    auxdata.MA(i).Joint(:,:) = DatStore.dM(:,i,:);  % moment arms
end
auxdata.DOFNames = DatStore.DOFNames;   % names of dofs

tau_act = 0.015; auxdata.tauAct = tau_act * ones(1,auxdata.NMuscles);       % activation time constant (activation dynamics)
tau_deact = 0.06; auxdata.tauDeact = tau_deact * ones(1,auxdata.NMuscles);  % deactivation time constant (activation dynamics)
auxdata.b = 0.1;                                                            % parameter determining transition smoothness (activation dynamics)

% Parameters of active muscle force-velocity characteristic
load ActiveFVParameters.mat
Fvparam(1) = 1.475*ActiveFVParameters(1); Fvparam(2) = 0.25*ActiveFVParameters(2);
Fvparam(3) = ActiveFVParameters(3) + 0.75; Fvparam(4) = ActiveFVParameters(4) - 0.027;
auxdata.Fvparam = Fvparam;

% Parameters of active muscle force-length characteristic
load Faparam.mat
auxdata.Faparam = Faparam;

% Parameters of passive muscle force-length characteristic
e0 = 0.6; kpe = 4; t50 = exp(kpe * (0.2 - 0.10e1) / e0);
pp1 = (t50 - 0.10e1); t7 = exp(kpe); pp2 = (t7 - 0.10e1);
auxdata.Fpparam = [pp1;pp2];
auxdata.ATendon=Misc.ATendon;

% Problem bounds
a_min = 0.01; a_max = 1;               % bounds on muscle activation
vA_min = -1/100; vA_max = 1/100;    % bounds on derivative of muscle activation (scaled)
F_min = 0; F_max = 5;               % bounds on normalized tendon force
dF_min = -100; dF_max = 100;        % bounds on derivative of normalized tendon force

% Time bounds
t0 = DatStore.time(1); tf = DatStore.time(end);
bounds.phase.initialtime.lower = t0; bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; bounds.phase.finaltime.upper = tf;

% Controls bounds
vAmin = vA_min./auxdata.tauDeact; vAmax = vA_max./auxdata.tauAct;
dFMin = dF_min*ones(1,auxdata.NMuscles); dFMax = dF_max*ones(1,auxdata.NMuscles);
aTmin = -1*ones(1,auxdata.Ndof); aTmax = 1*ones(1,auxdata.Ndof);
ControlExoMax = S.MaxTMult ;  ControlExoMin = S.MinTMult ;
bounds.phase.control.lower = [vAmin aTmin dFMin ControlExoMin]; bounds.phase.control.upper = [vAmax aTmax dFMax ControlExoMax];

% States bounds
actMin = a_min*ones(1,auxdata.NMuscles); actMax = a_max*ones(1,auxdata.NMuscles);
F0min = F_min*ones(1,auxdata.NMuscles); F0max = F_max*ones(1,auxdata.NMuscles);
Ffmin = F_min*ones(1,auxdata.NMuscles); Ffmax = F_max*ones(1,auxdata.NMuscles);
FMin = F_min*ones(1,auxdata.NMuscles); FMax = F_max*ones(1,auxdata.NMuscles);
bounds.phase.initialstate.lower = [actMin, F0min]; bounds.phase.initialstate.upper = [actMax, F0max];
bounds.phase.state.lower = [actMin, FMin]; bounds.phase.state.upper = [actMax, FMax];
bounds.phase.finalstate.lower = [actMin, Ffmin]; bounds.phase.finalstate.upper = [actMax, Ffmax];
% Integral bounds
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 100000000*(tf-t0);

% Path constraints
HillEquil = zeros(1, auxdata.NMuscles);
ID_bounds = zeros(1, auxdata.Ndof);
act1_lower = zeros(1, auxdata.NMuscles);
act1_upper = inf * ones(1, auxdata.NMuscles);
act2_lower = -inf * ones(1, auxdata.NMuscles);
act2_upper =  ones(1, auxdata.NMuscles)./auxdata.tauAct;
ExoPower_lower = -S.MaxPower;   ExoPower_upper = S.MaxPower;
bounds.phase.path.lower = [ID_bounds,HillEquil,act1_lower,act2_lower ExoPower_lower]; bounds.phase.path.upper = [ID_bounds,HillEquil,act1_upper,act2_upper ExoPower_upper];

% Eventgroup
% Impose mild periodicity
pera_lower = -1 * ones(1, auxdata.NMuscles); pera_upper = 1 * ones(1, auxdata.NMuscles);
perFtilde_lower = -1 * ones(1, auxdata.NMuscles); perFtilde_upper = 1 * ones(1, auxdata.NMuscles);
bounds.eventgroup.lower = [pera_lower perFtilde_lower]; bounds.eventgroup.upper = [pera_upper perFtilde_upper];

% Initial guess
N = length(DatStore.time);
guess.phase.time = DatStore.time;

% guess is currently not based on static optimzation => change this !
guess.phase.control = [zeros(N,auxdata.NMuscles) zeros(N,auxdata.Ndof) 0.01*ones(N,auxdata.NMuscles) zeros(N,auxdata.nExo)];
guess.phase.state =  [0.2*ones(N,auxdata.NMuscles) 0.2*ones(N,auxdata.NMuscles)];
guess.phase.integral = 0;

% Spline structures
for dof = 1:auxdata.Ndof
    for m = 1:auxdata.NMuscles
        auxdata.JointMASpline(dof).Muscle(m) = spline(DatStore.time,auxdata.MA(dof).Joint(:,m));
    end
    auxdata.JointIDSpline(dof) = spline(DatStore.time,DatStore.T_exp(:,dof));
    auxdata.JointIKSpline(dof) = spline(DatStore.time,DatStore.q_exp(:,dof).*pi./180);
end

for m = 1:auxdata.NMuscles
    auxdata.LMTSpline(m) = spline(DatStore.time,DatStore.LMT(:,m));
end


% GPOPS setup
setup.name = 'OptimalAnkleExo';
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'ma57';
setup.derivatives.derivativelevel = 'second';
setup.nlp.ipoptoptions.tolerance = 1e-6;
setup.nlp.ipoptoptions.maxiterations = 10000;
setup.derivatives.supplier = 'adigator';
setup.scales.method = 'none';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-3;
setup.mesh.maxiterations = 0;
setup.mesh.colpointsmin = 3;
setup.mesh.colpointsmax = 10;
setup.method = 'RPM-integration';
setup.displaylevel = 2;
NMeshIntervals = round((tf-t0)*Misc.Mesh_Frequency);
setup.mesh.phase.colpoints = 3*ones(1,NMeshIntervals);
setup.mesh.phase.fraction = (1/(NMeshIntervals))*ones(1,NMeshIntervals);
setup.functions.continuous = @Cont_nExo;
setup.functions.endpoint = @musdynEndpoint_FtildeState;

% ADiGator setup
input.auxdata = auxdata;
tdummy = guess.phase.time;
splinestruct = SplineInputData(tdummy,input);
splinenames = fieldnames(splinestruct);
for Scount = 1:length(splinenames)
    secdim = size(splinestruct.(splinenames{Scount}),2);
    splinestructad.(splinenames{Scount}) = adigatorCreateAuxInput([Inf,secdim]);
    splinestruct.(splinenames{Scount}) = zeros(0,secdim);
end
setup.auxdata.splinestruct = splinestructad;
if S.CreateAdigatorFiles
    adigatorGenFiles4gpops2(setup)
else
    clear Cont_nExoADiGatorGrd Cont_nExoADiGatorHes Cont_nExo
end

setup.functions.continuous =    @Wrap4_Cont_nExo;
setup.adigatorgrd.continuous =  @Cont_nExoGrdWrap;
setup.adigatorgrd.endpoint   =  @musdynEndpoint_FtildeStateADiGatorGrd;
setup.adigatorhes.continuous =  @Cont_nExoHesWrap;
setup.adigatorhes.endpoint   =  @musdynEndpoint_FtildeStateADiGatorHes;

end

