clear all;close all;clc

%%
[SNR0,h(1),RHO0] = cal_SNR2(selectedpair([1  7  13 14 ],1));
[SNR1,h(2),RHO1] = cal_SNR2(selectedpair([1  7  13 14 ],2));

% saveas(h(1),'day0sig.fig')
% saveas(h(2),'day1sig.fig')
% close all;
%%
hh1 = figure('Position',[1000 923 1000 400]);
hold on
errorbar([1:6]-.2,mean(RHO0),std(RHO0)/sqrt(numel(RHO0)),...
    'k.','Linewidth',1,'CapSize',15);
errorbar([1:6]+.2,mean(RHO1),std(RHO1)/sqrt(numel(RHO0)),...
    'k.','Linewidth',1,'CapSize',15);
b1(1) = bar([1:3]-.2,mean(RHO0(:,1:3)),.4,'FaceColor',[1 0.5 0.5]);
b1(2) = bar([4:6]-.2,mean(RHO0(:,4:6)),.4,'FaceColor',[0.5 0.5 .5]);
b2(1) = bar([1:3]+.2,mean(RHO1(:,1:3)),.4,'FaceColor',[1 0.5 0.5]/2);
b2(2) = bar([4:6]+.2,mean(RHO1(:,4:6)),.4,'FaceColor',[0.5 0.5 .5]/2);
title('pairwise correlation')

legend([b1(1:2) b2(1:2)],{'Pre-','Pre-','Post-','Post-'},'Orientation','horizontal')
legend('boxoff')
set(gca,'XTick',[1:6])
set(gca,'XTickLabel',{'1s','2s','3s','1s','2s','3s'});
set(gca,'YTick',0:.05:.2)
ylim([0 .075])
%% running-correlation
hh2 = histplt3(SNR0.corrs(selectedpair(:,1)),SNR1.corrs(selectedpair(:,2)),[1 0.5 0.5])
title(hh2(2),'running correlation-during stim')
hh3 = histplt3(SNR0.corrn(selectedpair(:,1)),SNR1.corrn(selectedpair(:,2)))
title(hh3(2),'running correlation-during blank')
% hh2(3) = histplt3(SNR0.corr(selectedpair(:,1)),SNR1.corr(selectedpair(:,2)))
%%
saveas(hh1,'pairwise-corr.fig')
saveas(hh2,'run corr-stim.fig')
saveas(hh3,'run corr-noise.fig')
%%
hh4 = histplt3(SNR0.skew(selectedpair(:,1)),SNR1.skew(selectedpair(:,2)))
title(hh4(2),'Skewness')
saveas(gca,'skewness.fig')