% cd('\\mps-zfs\data1\jsun\2pdata\Enhancement\completed\084\ROI2\')

if ~exist('day1peakf')
[day1peakf,path1]=uigetfile('peakSI.mat','select peakSI matfile for day1eye');
[day2peakf,path2]=uigetfile('peakSI.mat','select peakSI matfile for day2eye');
[prefilename,path]=uigetfile('selected_*pairs.mat','Pick the pre-matched pairs');
end

day1peak = matfile(fullfile(path1,day1peakf));
day2peak = matfile(fullfile(path2,day2peakf));
load(fullfile(path,prefilename),'selectedpair','peak_selected')
%% Load sig file, velocity, and stimulus info

[day1sig,p1]=uigetfile('.signals','select signal file for day1eye');
[day2sig,p2]=uigetfile('.signals','select signal file for day2eye');

day1file = load(fullfile(p1,day1sig),'-mat');
day2file = load(fullfile(p2,day2sig),'-mat');

pos_ = strfind(day1sig,'_');
day1run = load(fullfile(p1,[day1sig(1:pos_(3)) 'ball.mat']),'velocity');
day2run = load(fullfile(p2,[day2sig(1:pos_(3)) 'ball.mat']),'velocity');
day1stim = load(fullfile(p1,day1sig(1:pos_(3)-1)),'info');
day2stim = load(fullfile(p2,day2sig(1:pos_(3)-1)),'info');

%% Plottings
plotrunrep(day1peak.matrix,day2peak.matrix);


sigplt_onepage(day1file.sig,day1run.velocity,day1stim.info.frame(1:2:end),day1stim.info.stimtype)
mmpeaks1 = prctile(day1file.sig,[5 99]);
[SNR1,Cor1] = sort((mmpeaks1(2,:)+1)./(mmpeaks1(1,:)+1),'descend');
selected1 = Cor1(SNR1>1.29);
sigplt_onepage(day1file.sigsp(:,selected1)*2,day1run.velocity,day1stim.info.frame(1:2:end),day1stim.info.stimtype,Cor1(selected1))

sigplt_onepage(day1file.sig(:,selected1),day1run.velocity,day1stim.info.frame(1:2:end),day1stim.info.stimtype,Cor1(selected1))
sigplt_onepage(day2file.sig(:,selected1),day2run.velocity,day2stim.info.frame(1:2:end),day2stim.info.stimtype,Cor1(selected1))
mmpeaks2 = prctile(day2file.sig,[5 99]);
[SNR2,Cor2] = sort((mmpeaks2(2,:)+1)./(mmpeaks2(1,:)+1),'descend');
selected2 = Cor2(SNR2>1.29);
sigplt_onepage(day2file.sigsp(:,selected2)*2,day2run.velocity,day2stim.info.frame(1:2:end),day2stim.info.stimtype,Cor2(selected2))
%,day2file.sigsp(:,selected4)


mmpeaks = prctile([day1file.sig;day2file.sig],[5 99]);
[SNR,Cor] = sort((mmpeaks(2,:)+1)./(mmpeaks(1,:)+1),'descend');
selected = Cor(SNR>=1.76);
sigplt_onepage([day1file.sig(:,selected);day2file.sig(:,selected)],[day1run.velocity day2run.velocity],...
    [day1stim.info.frame(1:2:end);day2stim.info.frame(1:2:end)+13800],[day1stim.info.stimtype; day2stim.info.stimtype],selected)

sigplt_onepage([day1file.sigsp(:,selected);day2file.sigsp(:,selected)],[day1run.velocity day2run.velocity],...
    [day1stim.info.frame(1:2:end);day2stim.info.frame(1:2:end)+13800],[day1stim.info.stimtype; day2stim.info.stimtype],selected)


sigplt_onepage(day1file.sig(:,selected),day1run.velocity,day1stim.info.frame(1:2:end),day1stim.info.stimtype,selected)
sigplt_onepage(day1file.sigsp(:,selected),day1run.velocity,day1stim.info.frame(1:2:end),day1stim.info.stimtype,selected)

sigplt_onepage(day2file.sig(:,selected),day2run.velocity,day2stim.info.frame(1:2:end),day2stim.info.stimtype,selected)
sigplt_onepage(day2file.sigsp(:,selected)*3,day2run.velocity,day2stim.info.frame(1:2:end),day2stim.info.stimtype,selected)
%%
both = intersect(selected1,selected2);
unique1 = setdiff(selected1,both);
unique2 = setdiff(selected2,both);

figure;imshow(fig1img/2)
hold on;
for jj=intersect(both,find(pair0))
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1 1],'LineWidth',2);
end

for jj=intersect(unique1,find(pair0))
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 0 0],'LineWidth',1);
end

for jj=intersect(unique2,find(pair0))
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[0 1 0],'LineWidth',1);
end
%%

n=size(selectedpair,1);
rep = size(day1peak.sigF,2);
win_sig = day1peak.win_sig;

day0.sig = day1file.sig(:,selectedpair(:,1));
day1.sig = day2file.sig(:,selectedpair(:,2));

%%

    sigplt_onepage(cat(1,day0.sig,day1.sig), velocity,stim,Cor(:,example)); %

%%
day0.baseline = day1peak.sigF(1:win_sig(1)-1,:,:,selectedpair(:,1));
day1.baseline = day2peak.sigF(1:win_sig(1)-1,:,:,selectedpair(:,2));

%0-1s
day0.resp = day1peak.sigF(win_sig(1):win_sig(1)*2-2,:,:,selectedpair(:,1));
day1.resp = day2peak.sigF(win_sig(1):win_sig(1)*2-2,:,:,selectedpair(:,2));

%1-2s
day0.resp2 = day1peak.sigF(win_sig(1)*2-1:win_sig(1)*3-3,:,:,selectedpair(:,1));
day1.resp2 = day2peak.sigF(win_sig(1)*2-1:win_sig(1)*3-3,:,:,selectedpair(:,2));

day0.resp3 = day1peak.sigF(win_sig(1)*3-2:win_sig(1)*4-4,:,:,selectedpair(:,1));
day1.resp3 = day2peak.sigF(win_sig(1)*3-2:win_sig(1)*4-4,:,:,selectedpair(:,2));

day0.resp4 = day1peak.sigF(win_sig(1)*4-3:win_sig(1)*5-5,:,:,selectedpair(:,1));
day1.resp4 = day2peak.sigF(win_sig(1)*4-3:win_sig(1)*5-5,:,:,selectedpair(:,2));

%4-5s
day0.resp5 = day1peak.sigF(win_sig(1)*5-4:win_sig(1)*6-6,:,:,selectedpair(:,1));
day1.resp5 = day2peak.sigF(win_sig(1)*5-4:win_sig(1)*6-6,:,:,selectedpair(:,2));

for i=1:n
    day0.respshuffle(:,:,:,i) = day0.resp(:,randperm(rep),:,i);
    day1.respshuffle(:,:,:,i) = day1.resp(:,randperm(rep),:,i);
    
    day0.resp5shuffle(:,:,:,i) = day0.resp5(:,randperm(rep),:,i);
    day1.resp5shuffle(:,:,:,i) = day1.resp5(:,randperm(rep),:,i);
    
end

%%
day0corr = corr(day0.sig);
day1corr = corr(day1.sig);
neworder1 = dendrogramgraph(day1corr,day0corr,'Day1 correlation matrix','Day0 correlation matrix')
%%
day0baseline = reshape(day0.baseline,[],n);
day1baseline = reshape(day1.baseline,[],n);
day0basecor = corr(day0baseline);
day1basecor = corr(day1baseline);
daysbasecor = corr(day0baseline,day1baseline);

day0resp = reshape(day0.resp5,[],n);
day1resp = reshape(day1.resp5,[],n);
day0respshuffle = reshape(day0.respshuffle,[],n);
day1respshuffle = reshape(day1.respshuffle,[],n);

day0respcor = corr(day0resp);
day1respcor = corr(day1resp);
daysrespcor = corr(day0resp,day1resp);
day0respshcor = corr(day0resp,day0respshuffle);
day1respshcor = corr(day1resp,day1respshuffle);
daysrespshcor = corr(day0respshuffle,day1respshuffle);
%%
day0resp5 = reshape(day0.resp5,[],n);
day1resp5 = reshape(day1.resp5,[],n);
day0resp5shuffle = reshape(day0.resp5shuffle,[],n);
day1resp5shuffle = reshape(day1.resp5shuffle,[],n);
day0resp5shcor = corr(day0resp5,day0resp5shuffle);
day1resp5shcor = corr(day1resp5,day1resp5shuffle);
daysresp5shcor = corr(day0resp5shuffle,day1resp5shuffle);
%%
figure('Position',[200 0 1200 1200]);

% subplot(3,3,1)
% correlationgraph(day0basecor(neworder1,neworder1),'Day0 baseline correlation')
% subplot(3,3,2)
% correlationgraph(day1basecor(neworder1,neworder1),'Day1 baseline correlation')
% subplot(3,3,3)
% correlationgraph(daysbasecor(neworder1,neworder1),'Day0&1 baseline correlation')

% 
% subplot(3,3,4)
% correlationgraph(day0resp5cor(neworder1,neworder1),'Day0 resp correlation')
% subplot(3,3,5)
% correlationgraph(day1resp5cor(neworder1,neworder1),'Day1 resp correlation')
% subplot(3,3,6)
% correlationgraph(daysresp5cor(neworder1,neworder1),'Day0&1 resp correlation')
% 
subplot(3,3,7)
correlationgraph(day0resp5shcor(neworder1,neworder1),'Day0 4-5s shuffled resp correlation')
subplot(3,3,8)
correlationgraph(day1resp5shcor(neworder1,neworder1),'Day1 4-5s shuffled resp correlation')
subplot(3,3,9)
correlationgraph(daysresp5shcor(neworder1,neworder1),'Day0&1 4-5s shuffled resp correlation')


%%
[~,pref_ori] = nanmax(nanmax(peak_selected(3:4,:,:)));
pref_ori = squeeze(pref_ori);
pref_ori(pref_ori ==17) = NaN;
pref_ori(pref_ori>8) = pref_ori(pref_ori>8)-8;

M=pref_ori(neworder1);
M(M>4) = M(M>4)-8;

distancemap = abs(M*ones(1,n)- (M*ones(1,n))');
distancemap(distancemap>4) = distancemap(distancemap>4)-8;
distancemap = abs(distancemap);
figure;
imagesc(distancemap)
colormap('jet'); % set the colorscheme
colorbar

figure;imagesc(M');colormap('jet')
%%
% [~,pref_ori] = nanmax(nanmax(peak_selected(3:4,1:16,:)));
% pref_ori = squeeze(pref_ori);
% for ii=1:n
%     amp(ii) = peak_selected(4,pref_ori(ii),ii)./peak_selected(2,pref_ori(ii),ii);
% end
% for ii=1:n
%     MI(ii) = 1- peak_selected(1,pref_ori(ii),ii)./peak_selected(2,pref_ori(ii),ii);
% end
% MI = %log2(nansum(peak_selected(2,:,:))./nansum(peak_selected(1,:,:)));

% 
amp = nanmean(peak_selected(4,:,:))./nanmean(peak_selected(2,:,:));
amp = squeeze(amp);
MI = nanmean(peak_selected(1,:,:))./nanmean(peak_selected(2,:,:));
MI = 1-squeeze(MI);

% f=fit((1:n)',log2(amp(neworder1)),'exp1','StartPoint',[2,-1]);
% plot(f,1:n,log2(amp(neworder1)));
figure;hold on;
plot(1:n,MI(neworder1),'*');
plot(1:n,amp(neworder1),'o')
title('Enhancement')
ylim([0 2])
legend off
%%
M= MI(neworder1);
distancemap = abs(M*ones(1,n)- (M*ones(1,n))');
figure;
imagesc(distancemap)
colormap('jet'); % set the colorscheme
colorbar

figure;imagesc(M');colormap('jet')


%%
neworder0 = dendrogramgraph(day0corr,'Day0correlation matrix')
correlationgraph(day1corr(neworder1,neworder1),'Day1 correlation matrix')


