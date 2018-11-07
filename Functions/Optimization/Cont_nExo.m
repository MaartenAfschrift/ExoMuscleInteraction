function phaseout = Cont_nExo(input)

% Get input data
NMuscles        = input.auxdata.NMuscles;
Ndof            = input.auxdata.Ndof;
tauAct          = input.auxdata.tauAct;
tauDeact        = input.auxdata.tauDeact;
params          = input.auxdata.params;
splinestruct    = input.auxdata.splinestruct;
numColPoints    = size(input.phase.state,1);
Texo            = input.auxdata.Texo;
nExo            = input.auxdata.nExo;

% Get controls
vA   = 100*input.phase.control(:,1:NMuscles);
aT  = input.phase.control(:,NMuscles+1:NMuscles+Ndof);
dFtilde  = 10*input.phase.control(:,NMuscles+Ndof+1:NMuscles*2+Ndof);
uExo = input.phase.control(:,NMuscles*2+Ndof+1:NMuscles*2+Ndof+nExo);

% Get states
a       = input.phase.state(:,1:NMuscles);
Ftilde = input.phase.state(:,NMuscles+1:end);

% PATH CONSTRAINTS
% Activation dynamics - De Groote et al. (2009)
act1 = vA + a./(ones(size(a,1),1)*tauDeact);
act2 = vA + a./(ones(size(a,1),1)*tauAct);

% Hill-equilibrium constraint
% [Hilldiff,F] = ForceEquilibrium_FtildeState(a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam,input.auxdata.Atendon);
% [Hilldiff,F,lM,lMo,lMtilde,FMltilde,vM,vMtilde,FMvtilde,Fce,Fpe,FM] = ForceEquilibrium_FtildeState_FullInfo(...
%     a,Ftilde,dFtilde,splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam,input.auxdata.Atendon);

tension = ones(numColPoints,1)*input.auxdata.tension;
ATendon = input.auxdata.ATendon;
[Hilldiff,F, Fce, Fiso, vMmax, massM, vM , lMo,FMltilde,lMtilde] = ForceEquilibrium_FtildeState_all_LinTS(a,Ftilde,dFtilde,...
    splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam,tension,ATendon);
% [Hilldiff,F, Fce, Fiso, vMmax, massM, vM , lMo,FMltilde,lMtilde] = ForceEquilibrium_FtildeState_all(a,Ftilde,dFtilde,...
%     splinestruct.LMT,splinestruct.VMT,params,input.auxdata.Fvparam,input.auxdata.Fpparam,input.auxdata.Faparam,tension);

% If % ID torque is used to control exo



% Moments constraint
Topt = input.auxdata.Topt_res;
Tdiff = zeros(numColPoints,Ndof);
PowerExo=zeros(numColPoints,nExo);
for dof = 1:Ndof
    T_exp=splinestruct.ID(:,dof);
    if any(dof == input.auxdata.IndExo)
        iExo = dof == input.auxdata.IndExo;
        T_Exo_IDc=T_exp.*input.auxdata.ID_rel(iExo);
        Te = Texo(iExo).*uExo(:,iExo) + T_Exo_IDc;
        PowerExo(:,iExo) = Te .* splinestruct.IK_dot(:,dof);
        T_exp = T_exp - Te;
    end
    index_sel=(dof-1)*(NMuscles)+1:(dof-1)*(NMuscles)+NMuscles;
    T_sim=sum(F.*splinestruct.MA(:,index_sel),2) + Topt*aT(:,dof);
    Tdiff(:,dof) =  (T_exp-T_sim);
end

% constraint on exoskeleton power
phaseout.path = [Tdiff Hilldiff act1 act2 PowerExo];

% DYNAMIC CONSTRAINTS
% Activation dynamics is implicit
% Contraction dynamics is implicit
phaseout.dynamics = [vA dFtilde];


% OBJECTIVE FUNCTION
w0 = input.auxdata.w0;      % activation squared
w1 = input.auxdata.w1;      % reserve actuators
w2 = input.auxdata.w2;      % vA
w3 = input.auxdata.w3;      % metabolic energy

% compute metabolic energy
% compute metabolic energy
exc         = a;
act         = a;
vMtilde_E   = vM./lMo;
musclemass  = massM;                
pctst       = ones(numColPoints,1)*input.auxdata.ST_ratio;    % update this metric
vcemax      = 10.*lMo;
Fiso        = FMltilde;
b           = input.auxdata.b_Metab;

% energy_total=zeros(size(a));
% for m=1:NMuscles
%     [energy_total(:,m)] = ...
%     getMetabolicEnergySmooth2016all(exc(:,m),act(:,m),lMtilde(:,m),vMtilde_E(:,m),vM(:,m),Fce(:,m), ...
%                                     musclemass(:,m),pctst(:,m),vcemax(:,m),Fiso(:,m), 75,b);
% end
% if w3>0
% [energy_total] = ...
%     getMetabolicEnergySmooth2016all_vect_v2(exc,act,lMtilde,vMtilde_E,vM,Fce, ...
%     musclemass,pctst,vcemax,Fiso, b);


[energy_total] = ...
    getMetabolicEnergySmooth2016all_vect(exc,act,lMtilde,vMtilde_E,vM,Fce, ...
    musclemass,pctst,vcemax,Fiso, b);

% add power of reserve actuators to energy function
ReservePower =  (Topt*aT.*splinestruct.IK_dot);

phaseout.integrand = w0.*sum(a.^2,2)+ w1.*sum(aT.^2,2)+ w2*sum((vA/100).^2,2)+ w3*sum(energy_total.^2,2) + w3*sum(ReservePower.^2,2)*4;
    
% else
%     phaseout.integrand = w0.*sum(a.^2,2)+ w1.*sum(aT.^2,2)+ w2*sum((vA/100).^2,2);
end







