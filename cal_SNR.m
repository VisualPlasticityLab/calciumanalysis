function [SNR,h] = cal_SNR
% calculate Signal to noise ratio and plot with velocity and stim

[eyefile,p]=uigetfile('.signals','select signal file for eye1&eye2');
load(fullfile(p,eyefile),'sig','sig1','sigsp','stim','velocity','stimtype','-mat');
if ~exist('stim','var')
    velocityALL = [];
    stimtypeALL = [];
    stimALL = [];
    nframe = 0;
    eyefns = splitfile(eyefile);
    for n = 1:numel(eyefns)
        eyefn = eyefns{n};
        pos_ = strfind(eyefn,'_');
        load(fullfile(p,[eyefn(1:pos_(3)) 'ball.mat']),'velocity');
        try
            velocityALL = [velocityALL velocity];
        catch
            velocity = preprocess(fullfile(p,eyefn(1:pos_(3)-1)));
        end
        
        load(fullfile(p,eyefn(1:pos_(3)-1)),'info');
        stimtype = info.stimtype;
        stim = info.frame;
        stimtypeALL = [stimtypeALL; stimtype];
        stimALL = [stimALL; stim + sum(nframe)];
        save(fullfile(p,eyefn),'velocity','stim','stimtype','-append');
        nframe(n) = info.max_idx+1;
    end
    velocity = velocityALL;
    stimtype = stimtypeALL;
    stim = stimALL;
    save(fullfile(p,eyefile),'stim','velocity','stimtype','-append');
end

if ~exist('sig1')
    casignaltrend(fullfile(p,eyefile));
    load(fullfile(p,eyefile),'sig1','-mat');
end  

SNR.skew = skewness(sig);
mmpeaks = prctile(sig,[5 50 99]); 
SNR.ca = (mmpeaks(end,:)+1)./(mmpeaks(1,:)+1);
vel = downsample(velocity,round(numel(velocity)/size(sig,1)));
[SNR.corr,SNR.p] = sigcorr(sig,vel',1);
[SNR.corr1,SNR.p1] = sigcorr(sig1,vel',1);
hhh=figure;
subplot('Position',[0.1 0.1 .8 .6]); 
scatter(SNR.corr,SNR.corr1)
xlabel('signal-velocity corr');
ylabel('signal-velocity corr-corrected');
subplot('Position',[ 0.1 0.75 .8 .15]); 
histplt2(SNR.corr,SNR.corr1)
legend({'SV corr','after correction'},'orientation','hortizontal');
legend('boxoff')
saveas(hhh,'signal-velocity correction')
% mmpeaks_sp = prctile(sigsp,[5 50 99]); 
% SNR.sp = (mmpeaks_sp(end,:)+1)./(mmpeaks_sp(1,:)+1);

% 
% [SNRsorted,Cor] = sort(SNR.ca,'descend');
% selected = Cor(SNRsorted>=1.5);
[~,Cor] = sort(abs(SNR.corr),'descend');
selected = Cor;

% h(1)=figure;
% hold on;
% scatter(SNR.ca,SNR.sp,'ko')
% scatter(SNR.ca(selected),SNR.sp(selected),'ro')
% xlabel('Calcium signal noise ratio 99%-5%')
% ylabel('Spikie signal noise ratio 99.5% - 50%');
% h(2) = sigplt_onepage(sigsp(:,selected),velocityALL,stimALL,stimtypeALL,selected)
save(fullfile(p,eyefile),'SNR','-append');

h = sigplt_onepage(sig(:,selected),velocity,stim,stimtype,selected)
 saveas(h,'traces_ca.fig')
% saveas(h(3),'traces_sp.fig')

