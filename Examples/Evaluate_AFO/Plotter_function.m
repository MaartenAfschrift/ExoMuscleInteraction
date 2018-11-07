%% Plotter evaluate AFO conditions

clear all; close all; clc;
MainDir     = 'C:\Users\u0088756\Documents\PostDoc\Software\ExoSimulation\BenchMarkTest';
Datapath    = fullfile(MainDir,'Data','Gait92Model');
S.OutPath   = fullfile(Datapath,'EvaluateAFO_results');
load(fullfile(S.OutPath,'EvaluateAFO.mat'));

mInds =32:34;


t=res.time-res.time(1);
lw=2;
% plot ID torque
figure();
subplot(1,4,1:2)
plot(t,squeeze(TID_Vect(:,5,1)),'k','LineWidth',lw); hold on;
plot(t,squeeze(Texo_Vect(:,3,1)),'r','LineWidth',lw); hold on;   % optimal assitance
plot(t,squeeze(Texo_Vect(:,3,2)),'g','LineWidth',lw); hold on;   % ID assistance
plot(t,squeeze(Texo_Vect(:,3,3)),'b','LineWidth',lw); hold on;   % ID resistance
delete_box


% change metabolic energy
% figure(); 

TotPower=squeeze(sum(E_Vect,2)); 
MetabEnergy=zeros(1,4);
for i=1:4
	MetabEnergy(i) =  trapz(res.time,TotPower(:,i));
end

% energy change
E_change = -MetabEnergy+MetabEnergy(4);
E_change_rel=E_change./MetabEnergy(4);
subplot(1,4,3)
bar(1,E_change_rel(1),'r'); hold on;
bar(2,E_change_rel(2),'g');
bar(3,E_change_rel(3),'b');

% exo efficiency
ExoWork = PosWork_Vect(end,:)-PosWork_Vect(end,4);
efficiency = (E_change*0.25)./ExoWork;

subplot(1,4,4)
bar(1,efficiency(1),'r'); hold on;
bar(2,efficiency(2),'g');
bar(3,efficiency(3),'b');


% 
% figure();
% subplot(3,1,1)
% bar(E_change)
% subplot(3,1,2)
% bar(ExoWork)
% subplot(3,1,3)
% bar(efficiency)
% 
% figure();
% i=3;
% plot(res.time,squeeze(lMtilde_Vect(:,mInds(i),1)),'r'); hold on;
% plot(res.time,squeeze(lMtilde_Vect(:,mInds(i),2)),'g'); hold on;
% plot(res.time,squeeze(lMtilde_Vect(:,mInds(i),3)),'b'); hold on;
% plot(res.time,squeeze(lMtilde_Vect(:,mInds(i),4)),'k'); hold on;
% % end
% 
% figure();
% for i=1:3
%     subplot(1,3,i)
%     plot(res.time,squeeze(MusclePower_Vect(:,mInds(i),1)),'r'); hold on;
%     plot(res.time,squeeze(MusclePower_Vect(:,mInds(i),2)),'g'); hold on;
%     plot(res.time,squeeze(MusclePower_Vect(:,mInds(i),3)),'b'); hold on;
%     plot(res.time,squeeze(MusclePower_Vect(:,mInds(i),4)),'k'); hold on;
% end
