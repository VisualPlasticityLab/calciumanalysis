% upper = load('SST01_1031_001&SST01_1031_002_1&SST01_1107_101&SST01_1107_103_2_selected10')
% lower = load('SST01_1031_001&SST01_1031_002_2&SST01_1107_101&SST01_1107_103_3_selected18')


upper = load('2018_627_001&2018_627_002_1&2018_703_001&2018_703_002_1_selected18');
middle = load('2018_627_001&2018_627_002_2&2018_703_001&2018_703_002_2_selected32');
lower = load('2018_627_001&2018_627_002_3&2018_703_001&2018_703_002_3_selected29');

ODI1s = cat(2,upper.ODI1s,middle.ODI1s,lower.ODI1s);
ODI1r = cat(2,upper.ODI1r,middle.ODI1s,lower.ODI1r);
numcell = num2str(size(ODI1s,2));
filenm = ['combined_selected' numcell]

%%
sets = 2; %define how many columns
h(2)=figure('Position',[200 200 300*sets 500]);

subplot(2,sets,1);title('STILL','FontSize',12)
hold on
plot([-1 1],[-1 1],'--','Color',[0.5 0.5 .5 .5])
% scatter(ODI1s(1,:),ODI1s(2,:),peak*50,'ko','filled'); 
scatter(ODI1s(1,:),ODI1s(2,:),'ko','filled'); 
scatter(nanmean(ODI1s(1,:)),nanmean(ODI1s(2,:)),100,'kx');
alpha(.3)
xlabel(['ODIbefore \mu= ' sprintf('%.2f',nanmean(ODI1s(1,:)))])
ylabel(['ODIafter \mu= ' sprintf('%.2f',nanmean(ODI1s(2,:)))])
axis square

subplot(2,sets,sets+1);hold on;
histplt2(ODI1s(1,:),ODI1s(2,:),-1.0:2/11:1.0);
legend('before','after','Location','Northwest')
legend('Boxoff')
xlabel('ODI@best direction')
ylabel('Cell count')
title(['N=' numcell],'FontSize',12)
axis square

% RUN  condition
subplot(2,sets,2);title('RUN','FontSize',12)
hold on
plot([-1 1],[-1 1],'--','Color',[0.5 0.5 .5 .5])
% scatter(ODI1r(1,:),ODI1r(2,:),peak*50,'ko','filled'); 
scatter(ODI1r(1,:),ODI1r(2,:),'ko','filled'); 
scatter(nanmean(ODI1r(1,:)),nanmean(ODI1r(2,:)),100,'kx');
alpha(.3)

xlabel(['ODIbefore \mu= ' sprintf('%.2f',nanmean(ODI1r(1,:)))])
ylabel(['ODIafter \mu= ' sprintf('%.2f',nanmean(ODI1r(2,:)))])
axis square

subplot(2,sets,sets+2);hold on;
histplt2(ODI1r(1,:),ODI1r(2,:),-1.0:2/11:1.0);
legend('before','after','Location','Northwest')
legend('Boxoff')
xlabel('ODI@best direction')
ylabel('Cell count')
title(['N=' numcell],'FontSize',12)
axis square

%%
saveas(h(2),[filenm '.png'])
save(filenm,'ODI1r','ODI1s')
close all