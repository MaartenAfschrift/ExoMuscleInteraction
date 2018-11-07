%% Evaluate simple results Exo
%-----------------------------

% path information
MainDir     = 'C:\Users\u0088756\Documents\PostDoc\Software\ExoSimulation\BenchMarkTest';
Datapath    = fullfile(MainDir,'Data','Gait92Model');
load(fullfile(Datapath,'results','Energy_TendonSim.mat'));


nSim=9;

mInds =32:34;
tv = linspace(0,100,length(res.time));
GrayCol=hsv(nSim);

%% plot exo torque
h=figure();
for i=1:nSim
    subplot(3,3,1)
    plot(tv,squeeze(Texo_Vect(:,3,i)),'Color',GrayCol(i,:)); hold on;
end
subplot(3,3,1)
plot(tv,squeeze(TID_Vect(:,5,i,1)),'--k','LineWidth',2); hold on;
xlabel('% gait cycle');     ylabel('Ankle torque [Nm]');
set(gcf,'Name','Complex-Actuator');

%% Plot muscle activity gastrocn and soleus
for i=1:nSim
    subplot(3,3,2)
    plot(tv,squeeze(MActivation_Vect(:,mInds(1),i)),'Color',GrayCol(i,:)); hold on;
    subplot(3,3,3)
    plot(tv,squeeze(MActivation_Vect(:,mInds(3),i)),'Color',GrayCol(i,:)); hold on;
end
for i=1:2
    subplot(3,3,1+i)
    xlabel('% gait cycle'); 
    ylabel('Muscle Activity []');
end

%% Plot fiber length

for i=1:nSim
    subplot(3,3,4)
    plot(tv,squeeze(lMtilde_Vect(:,mInds(1),i)),'Color',GrayCol(i,:)); hold on;
    plot(tv,squeeze(lMtilde_Vect2(:,mInds(1),i)),'--k'); hold on;
    subplot(3,3,5)
    plot(tv,squeeze(lMtilde_Vect(:,mInds(3),i)),'Color',GrayCol(i,:)); hold on;
    plot(tv,squeeze(lMtilde_Vect2(:,mInds(3),i)),'--k'); hold on;
end
for i=1:2
    subplot(3,3,3+i)
    xlabel('% gait cycle'); 
    ylabel('Norm fiber length []');
end

for i=1:nSim
    subplot(3,3,6)
    plot(tv(1:end-3),squeeze(MusclePower_Vect(1:end-3,mInds(1),i)),'Color',GrayCol(i,:)); hold on;
%     plot(tv(1:end-3),squeeze(MusclePower_Vect2(1:end-3,mInds(1),i)),'--k'); hold on;
    subplot(3,3,7)
    plot(tv(1:end-3),squeeze(MusclePower_Vect(1:end-3,mInds(3),i)),'Color',GrayCol(i,:)); hold on;
%     plot(tv(1:end-3),squeeze(MusclePower_Vect2(1:end-3,mInds(3),i)),'--k'); hold on;
end
for i=1:2
    subplot(3,3,5+i)
    xlabel('% gait cycle'); 
    ylabel('Muscle Power (W)');
end

%% Plot energy consumption
for i=1:nSim
    subplot(3,3,8)
    E=sum(E_Vect(:,:,i),2); Etot = trapz(res.time,E);
    plot(StiffnessVect(i),Etot,'*','Color',GrayCol(i,:)); hold on;
     plot(StiffnessVect(i),Etot,'O','Color',GrayCol(i,:)); hold on;
    E=sum(E_Vect2(:,:,i),2); Etot = trapz(res.time,E);
    plot(StiffnessVect(i),Etot,'Ok'); hold on;
   
end

subplot(3,3,8)
xlabel('Stiffness');     ylabel('Metab energy [J]');
% set(gca,'YLim',[100 400]);
delete_box

%% Plot percentage change in Energy
P_Vect=squeeze(sum(E_Vect,2));
P_Vect2=squeeze(sum(E_Vect2,2));
subplot(3,3,9);
for i=1:nSim
    Wp=trapz(res.time,P_Vect(:,i));
    Wunp=trapz(res.time,P_Vect2(:,i));
    Gain=(Wunp-Wp)./Wunp;
    b=bar(i,Gain);
    set(b,'FaceColor',GrayCol(i,:)); hold on;
end
ylabel('Perc. Reduction Metab. E');

