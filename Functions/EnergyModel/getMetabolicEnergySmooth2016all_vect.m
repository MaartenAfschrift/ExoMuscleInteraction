% This function computes the muscle energy expenditure based on the model
% of Uchida2016. Energy_total is in W. 
% This function provides a smooth approximation of the original function
% for use in optimization problems.
% As compared to getMetabolicEnergySmoothUmberger2010all, this versions 1)
% includes negative mechanical work, 2) ensures that the total power is 
% non-negative, 3) use an order recruitment model
%
% Code inspired from Maarten Afschrift and Azin Zargham
% Date: 14/02/2018
% Inputs:
%   act: muscle activation
%   lMtilde: normalized muscle fiber length
%   vMtilde: normalized muscle fiber velocity (positive is lengthening).
%   Note: defined for this model as vM/lMopt, which is a different metric
%   to the typical vM/vMax
%   vM: muscle fiber velocity (positive is lengthening)
%   Fce: muscle force from contractile element (length + velocity
%   components but not passive component)
%   musclemass: mass of muscle, mass = PCSA*rho*lMopt, PCSA = Fmax/sigma
%   where rho is density (1058.7 kg/m³) and sigma is tension (0.25MPa)
%   pctst: percentage of slow twitch fibers (0-1)
%   vcemax: maximal muscle fiber velocity (default is 10*lMopt)
%   Fiso: normalized muscle force from active force-length relationship
%   modelmass: mass model
%   b: parameter determining transition smoothness for tanh approximations
% Outputs:
%   energy_total: totat metabolic energy (W)
%   energy_am: energy from activation and maintenance (W)
%   energy_sl: energy from shortening and lengthening (W)
%   energy_mech: mechanical energy (W)
%   energy_model: energy model (includes basal rate) (W)
% Note: the energy_total might be different than the sum of the different
% component if the total heat rate was clamped to one.

function [energy_total,energy_am,energy_sl,energy_mech] = ...
    getMetabolicEnergySmooth2016all_vect(exc,act,lMtilde,vMtilde,vM,Fce, ...
                                    musclemass,pctst,vcemax,Fiso, b)
%% Order recruitment model
% The percentage of slow-twitch fibers recruited is function of excitation
% f_rec_slow = getPctSTSmooth_vect(exc,pctst,b);      
f_rec_slow = pctst;
                                
%% Some parameters
% Percentage of slow twitch (st) and fast twitch (ft) fibers
f_rec_fast = 1-f_rec_slow;
% Scaling factor: aerobic (1.5) or anaerobic (1) activities
s = 1.5;            

%% Excitations and activations
% a = exc + (-exc+act)/2.*(0.5+0.5*tanh(b*(act-exc)));
a = exc;
%% Heat activation and maintenanance
% lMtilde is on the descending leg of the active force-length relationship
% i.e. lMtilde > 1 || lM > lMopt
lMtilde_des = 0.5 + 0.5*tanh(b*(lMtilde-1.*ones(size(lMtilde))));        
hdotam = 128*f_rec_fast+25 + (-0.6*(128*f_rec_fast+25) ...
         + 0.6*(128*f_rec_fast+25).*Fiso).*lMtilde_des;       

aam         = a.^(0.6);
energy_am   = hdotam.*aam*s;

%% Heat shortening and lengthening
% vM is positive (muscle lengthening)
vMtilde_pos = 0.5 + 0.5*tanh(b*(vMtilde));
% vM is negative (muscle shortening)
vMtilde_neg = 1-vMtilde_pos;

coef_hs_ft  = 1*153.*ones(size(vcemax))./vcemax;
coef_hs_st  = 4*25.*ones(size(vcemax))./(vcemax/2.5); 
coef_hl     = 4*coef_hs_st; % different compared to Umberger2010

hdotsl = coef_hl.*vMtilde + (-coef_hl.*vMtilde ...
    - coef_hs_st.*vMtilde.*f_rec_slow ...
    - coef_hs_ft.*vMtilde.*f_rec_fast).*vMtilde_neg;

energy_sl = (hdotsl.*a.*s + ...
    (-hdotsl.*a.*s + hdotsl.*(a.^2)*s).*vMtilde_neg).*(1-lMtilde_des) + ...
    (hdotsl.*a.*s.*Fiso + (-hdotsl.*a.*s.*Fiso + ...
    hdotsl.*(a.^2)*s.*Fiso).*vMtilde_neg).*(lMtilde_des);

%% Mechanical work
% Negative mechanical work allowed
energy_mech = -1*Fce.*vM; % different compared to Umberger2010

%% Total power is non-negative
% If necessary, increase the shortnening heat rate
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Umberger2010MuscleMetabolicsProbe.cpp#L410
Edot_Wkg_beforeClamp = energy_am + energy_sl + energy_mech;
Edot_Wkg_beforeClamp_neg = 0.5+(0.5*tanh(b*(-Edot_Wkg_beforeClamp)));
energy_sl = energy_sl - Edot_Wkg_beforeClamp.*Edot_Wkg_beforeClamp_neg;

%% The total rate of energy liberation is not allowed to drop below 1(W/kg)
% https://github.com/opensim-org/opensim-core/blob/master/OpenSim/Simulation/Model/Umberger2010MuscleMetabolicsProbe.cpp#L432
totalHeatRate = energy_am + energy_sl;
totalHeatRate = totalHeatRate + (-totalHeatRate + ...
    ones(size(totalHeatRate))).*(0.5 + 0.5*tanh(b*(1-totalHeatRate)));

%% Account for muscle mass
energy_am = energy_am.*musclemass;
energy_sl = energy_sl.*musclemass;
energy_total = totalHeatRate.*musclemass+energy_mech;


end

