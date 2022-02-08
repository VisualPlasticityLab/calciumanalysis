function corrplt(pre,post,thr)
run_corr = corr(pre.velocity',post.velocity');
% nonrun=find(pre.velocity<=thr&post.velocity<=thr);
section=(pre.velocity<=thr&post.velocity<=thr);
sect = (diff(section)==1);
% nonrun_sect = 
for n=1:npair
    run1_corr(n,:) = corr(pre.spks(n,section)',pre.velocity(section)');
    run2_corr(n,:) = corr(post.spks(n,section)',post.velocity(section)');
    sessions_corr(n,:) = corr(pre.spks(n,section)',post.spks(n,section)');
end
% Plot correlation
figure;
boxplot([run1_corr, sessions_corr, run2_corr]);
hold on;
plot([1;2;3], [run1_corr, sessions_corr,run2_corr],'x-','Color',[.5 .5 .5 .3],'LineWidth',.5);
%plot(1.5,sessions_corr)
xlim([.5 3.5])
xticks(1:3)
xticklabels({'run-pre','pre-post','run-post'})
title(sprintf('correlation-allpairs,runcorr=%.2f',run_corr))

