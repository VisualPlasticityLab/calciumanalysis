%%%% load SI. mat and plot responses


for ori=1:9;
figure;hold on;title([num2str(45*(ori-1)) 'degree']);axis equal
scatter(squeeze(SI.peakS(:,ori,:)),squeeze(SI.peakR(:,ori,:)),'o');%,'jitter','on', 'jitterAmount', 0.10);
plot([0.01 1],[0.01 1],'r--');
xlabel('still');ylabel('running');
% h=gca; h.XScale='log'; h.YScale='log';
end

figure;hold on;title('Running-Still');axis equal
MI2=(squeeze(SI.peakR(:,2,:))-squeeze(SI.peakS(:,2,:)));
MI6=(squeeze(SI.peakR(:,6,:))-squeeze(SI.peakS(:,6,:)));
scatter(MI2,MI6,'o');%,'jitter','on', 'jitterAmount', 0.10);
plot([-1 1 ],[-1 1],'r--');

figure;
runstim_R=SI.peakR(:,3,:)+SI.peakR(:,11,:);
runstim_S=SI.peakS(:,3,:)+SI.peakS(:,11,:);
stillstim_R= SI.peakR(:,7,:)+SI.peakR(:,15,:);
stillstim_S= SI.peakS(:,7,:)+SI.peakS(:,15,:);
plot(squeeze(cat(2,runstim_R,stillstim_R,runstim_S,stillstim_S)),'o-');
xlim([0.5 4.5])
g=gca;
g.XTick=1:4;
g.XTickLabel={'runstim_R','stillstim_R','runstim_S','stillstim_S'}
ylabel('Response Amp dF/F')
saveas(gcf,'Drft and State')
saveas(gcf,'Drft and State.png')

figure;
runstim_R=SI.peakR(:,3,:)+SI.peakR(:,11,:);
stillstim_R= SI.peakR(:,7,:)+SI.peakR(:,15,:);
runstim_S=SI.peakS(:,3,:)+SI.peakS(:,11,:);
stillstim_S= SI.peakS(:,7,:)+SI.peakS(:,15,:);
plot(squeeze(cat(2,runstim_R+runstim_S,stillstim_R+stillstim_S)),'o-');
xlim([.5 2.5])
g=gca;
g.XTick=1:2;
g.XTickLabel={'runstim_R_S','stillstim_R_S'}
ylabel('Response Amp dF/F')
saveas(gcf,'Drft only')
saveas(gcf,'Drft only.png')