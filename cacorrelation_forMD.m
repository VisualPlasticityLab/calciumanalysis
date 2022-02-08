cd('\\mps-zfs\data1\mdadarla\CA2p\6284\150530\m1')

%% Load sig file, velocity, and stimulus info
load('Stim_data')

[day1sig,p1]=uigetfile('.signals','select signal file for day1eye');
[day2sig,p2]=uigetfile('.signals','select signal file for day2eye');

day1file = load(fullfile(p1,day1sig),'-mat');
day2file = load(fullfile(p2,day2sig),'-mat');


pos_ = strfind(day1sig,'_');
day1run = load(fullfile(p1,[day1sig(1:pos_(3)) 'ball.mat']),'velocity');
day2run = load(fullfile(p2,[day2sig(1:pos_(3)) 'ball.mat']),'velocity');

multiplane = str2num(day1sig(pos_(3)+1: pos_(4)-1));
nplane = 3;

v1 = downsample(day1run.velocity,nplane,multiplane);
v2 = downsample(day2run.velocity,nplane,multiplane);
ttlonfr1 = floor((day1file.ttlonfr1-multiplane)/nplane)+1;
ttlonfr2 = floor((day2file.ttlonfr2-multiplane)/nplane)+1;
Ampl1 = day1file.Ampl1;
Ampl2 = day2file.Ampl2;
sig1 = day1file.sig(:,1:end-1);
sig2 = day2file.sig(:,1:end-1);
%% Plottings
% plotrunrep(day1peak.matrix,day2peak.matrix);
% sigplt_onepage(sig1,day1run.velocity,ttlonfr1,Ampl1)
d1sk = skewness(sig1);
d2sk = skewness(sig2);
dsk = skewness([sig1;sig2]);

selected1 = d1sk>1.25;
selected2 = d2sk>1.25;
selected0 = (d1sk>1.25)|(d2sk>1.25);
selected = dsk>1.25;

Cor = 1:numel(d1sk);

level = 6;
C = discretize(dsk,level)*20;
figure('Position',[341 848 1219 490]);
subplot(1,2,1)
hold on
scatter(d1sk,d2sk,C)
plot([-2 10],[-2 10],'k--')
axis square
x = [-2 1.25 1.25 -2];
y = [-2 -2 1.25 1.25];
patch(x',y','k')
alpha(.3)
title('Signal skewness')
xlabel('before')
ylabel('after')
subplot(1,2,2);
histplt2(d1sk,d2sk)

figure;hold on;scatter(d1sk,d2sk-d1sk)
plot([-2 10],[0 0],'k--')
plot([1.25 1.25],[-6 8],'k--')

sigplt_onepage(sig1(:,selected),v1,ttlonfr1,Ampl1,Cor(selected))
sigplt_onepage(sig2(:,selected),v2,ttlonfr2,Ampl2,Cor(selected))

sigplt_onepage([sig1(:,selected);sig2(:,selected)],[v1 v2],...
    [ttlonfr1;ttlonfr2+23460/3],[Ampl1; Ampl2],Cor(selected))
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
day0.sig = sig1(:,selected0);
day1.sig = sig2(:,selected0);
n=size(day0.sig,2);
%%
day0.baseline = [];
day1.baseline = [];

for ii = 2:numel(ttlonfr1)
    day0.baseline(end+1:end+38,:) = day0.sig(ttlonfr1(ii)-40:ttlonfr1(ii)-3,:);
    day1.baseline(end+1:end+38,:) = day1.sig(ttlonfr2(ii)-40:ttlonfr2(ii)-3,:);
    
end

day0.resp = [];
day1.resp = [];
%0-1s
for  ii = 1:numel(ttlonfr1)
    day0.resp(end+1:end+11,:) = day0.sig(ttlonfr1(ii):ttlonfr1(ii)+10,:);
    day1.resp(end+1:end+11,:) = day1.sig(ttlonfr2(ii):ttlonfr2(ii)+10,:);
end

for i=1:n
    day0.respshuffle(:,:,:,i) = day0.resp(:,randperm(rep),:,i);
    day1.respshuffle(:,:,:,i) = day1.resp(:,randperm(rep),:,i);
end

%%
day0corr = corr(day0.sig);
day1corr = corr(day1.sig);
neworder1 = dendrogramgraph(day1corr,day0corr,'Day1 correlation matrix','Day0 correlation matrix')
%%
day0baseline = day0.baseline;
day1baseline = day1.baseline;
day0basecor = corr(day0baseline);
day1basecor = corr(day1baseline);
daysbasecor = corr(day0baseline,day1baseline);

day0resp = reshape(day0.resp,[],n);
day1resp = reshape(day1.resp,[],n);
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
nrow = 3;
ncol = 2;

subplot(ncol,3,1)
correlationgraph(day0basecor(neworder1,neworder1),'Day0 baseline correlation')
subplot(ncol,3,2)
correlationgraph(day1basecor(neworder1,neworder1),'Day1 baseline correlation')
subplot(ncol,3,3)
correlationgraph(daysbasecor(neworder1,neworder1),'Day0&1 baseline correlation')


subplot(ncol,3,4)
correlationgraph(day0respcor(neworder1,neworder1),'Day0 resp correlation')
subplot(ncol,3,5)
correlationgraph(day1respcor(neworder1,neworder1),'Day1 resp correlation')
subplot(ncol,3,6)
correlationgraph(daysrespcor(neworder1,neworder1),'Day0&1 resp correlation')
% % 
% subplot(3,3,7)
% correlationgraph(day0resp5shcor(neworder1,neworder1),'Day0 4-5s shuffled resp correlation')
% subplot(3,3,8)
% correlationgraph(day1resp5shcor(neworder1,neworder1),'Day1 4-5s shuffled resp correlation')
% subplot(3,3,9)
% correlationgraph(daysresp5shcor(neworder1,neworder1),'Day0&1 4-5s shuffled resp correlation')


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


