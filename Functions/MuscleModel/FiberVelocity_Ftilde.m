% This function computes the muscle fiber length from the normalized tendon
% force

function [lM,lMtilde,vM,vMtilde,lTtilde ] = FiberVelocity_Ftilde(Ftilde,dfse,params,lMT,vMT,Atendon,bool_LinearSpring)


% input arguments
lMo = ones(size(Ftilde,1),1)*params(2,:);
lTs = ones(size(Ftilde,1),1)*params(3,:);
alphao = ones(size(Ftilde,1),1)*params(4,:);
vMmax = ones(size(Ftilde,1),1)*params(5,:);

% Non-linear tendon
if bool_LinearSpring
    lTtilde = Ftilde./Atendon+1;
else
    warning('Non linear spring not implemented - using linear spring')
    lTtilde = Ftilde./Atendon+1;
end

% Hill-model relationship
lM = sqrt((lMo.*sin(alphao)).^2+(lMT-lTs.*lTtilde).^2);
lMtilde = lM./lMo;
if bool_LinearSpring
    vT = lTs.*dfse./Atendon;
else
    vT = lTs.*dfse./Atendon;
end
cos_alpha = (lMT-lTs.*lTtilde)./lM;
vM = (vMT-vT).*cos_alpha;
vMtilde = vM./vMmax;
end

