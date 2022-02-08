function [SNR,h,RHO]=cal_SNR2(selected)
%% load/save stim, velocity, stim etc info
[eyefile,p]=uigetfile('.signals','select signal file');
load(fullfile(p,eyefile),'sig','stim','velocity','stimtype','-mat');
if ~exist('stim')
    pos_ = strfind(eyefile,'_');
    load(fullfile(p,[eyefile(1:pos_(3)) 'ball.mat']),'velocity');
    load(fullfile(p,eyefile(1:pos_(3)-1)),'info');
    stimtype = info.stimtype;
    stim = info.frame;
    save(fullfile(p,eyefile),'velocity','stim','stimtype','-append');
end
%% Plot selected
h = sigplt_onepage(sig(:,selected),velocity,stim,stimtype,selected);
% saveas(h,'traces_ca.fig')
% 
% RHOall = corr(sig(:,selected));
% nanmean(RHOall(~triu(ones(numel(selected)))))
st = ceil(stim/3);
stON = st(1:2:end-1);
stOFF = st(2:2:end);
stONtime = min(stOFF-stON);

%% calculate pairwise correlation
% dur = 1:round(stONtime/3);
% for i=1:3
% resp = repmat(stON',numel(dur),1)+repmat(dur'+(i-1)*round(stONtime/3),1,numel(stON));
% noise = repmat(stOFF',numel(dur),1)+repmat(dur'+(i-1)*round(stONtime/3),1,numel(stON));
% RHOs = corr(sig(resp,selected));
% RHOn = corr(sig(noise,selected));
% RHO(:,i) = (RHOs(~triu(ones(numel(selected)))));
% RHO(:,i+3) = (RHOn(~triu(ones(numel(selected)))));
% end


dur = 1:round(stONtime/23);
resp = repmat(stON',numel(dur),1)+repmat(dur',1,numel(stON));
noise = repmat(stOFF',numel(dur),1)+repmat(dur',1,numel(stON));
RHO=[];
RHOs = corr(sig(resp,selected));
RHO(:,1) = (RHOs(~triu(ones(numel(selected)))));
RHOn = corr(sig(noise,selected));
RHO(:,2) = (RHOn(~triu(ones(numel(selected)))));
save(fullfile(p,eyefile),'RHO','-append');
%%  calculate running correlation
SNR.skew =skewness(sig);
vel = downsample(velocity,round(numel(velocity)/size(sig,1)));
[SNR.corr,SNR.p] = sigcorr(sig,vel',1);
title('sig velocity correlation-allcells')

SNR.corrs(:,1) = sigcorr(sig(resp,:),vel(1,resp)',1);
SNR.corrn(:,1) = sigcorr(sig(noise,:),vel(1,noise)',1);

save(fullfile(p,eyefile),'SNR','-append');



