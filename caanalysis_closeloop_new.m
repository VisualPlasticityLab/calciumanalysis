
%function ca_analysis;
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clearvars -except gd1;global info;
%% load calcium signal and open recorded image
[fname,p]=uigetfile('.signals','load calcium signal data .signals');
load(fullfile(p,fname),'-mat');

if ~exist('sig')
    sig=sig_chunk;
end
if size(sig,1)<size(sig,2)
    sig=sig';
end
ncell=size(sig,2);
%% sync with stimulus data and plot
try
    pos = strfind(strtok(fname,'c'),'_');
    fn = fname(1:pos(3)-1);
catch
    [temp,path]=uigetfile('*.sbx','Pick the corresponding sbx file');
    fname=fullfile(path,strtok(temp,'.'));
end
velocity=preprocess(fn); %velocity has ball info and info has stimulus info

%% 2. interpolate ca signal, derive smoothed or spiking version
Nframes = info.max_idx+1;
if info.volscan == 0 & Nframes/size(sig,1) < 2
%     x=1:size(sig,1);
        nframes = Nframes;
        nplane = 1;
else
%     nplane = info.otparam(3);
    multiplane = str2num(fname(pos(end)-1));
    nplane = round(Nframes/size(sig,1));
    x=multiplane:nplane:Nframes;
    xx=1:Nframes;
    % sig=interp1(x,full(sig),xx,'spline','extrap');
    frameorg = floor((info.frame-multiplane)/nplane)+1;
    velocity = interp1(xx,velocity,x,'previous');
    nframes = floor((Nframes-multiplane)/nplane)+1;
end
CAframeHz =info.resfreq/info.recordsPerBuffer/nplane;
CAlineHz = info.config.lines/nplane;
newtime=info.time(info.counting==1);
stimtypeorg =info.stimtype;
%just for plotting purpose,assign stimtype to all the ON frames
stim=zeros(1,nframes);
%for each stimtype, it sends signal to info.frame for its beginning and ending
for ii=1:numel(stimtypeorg)
    stim(frameorg(ii*2-1):frameorg(ii*2))=stimtypeorg(ii); 
end
%% CORRECTION! only count when the stimulus flipps
novel = find([true;stimtypeorg(2:end)~=stimtypeorg(1:end-1)]);
frame = zeros(numel(novel)*2,1);
frame(1:2:end) = frameorg(novel*2-1);
novelend = find([stimtypeorg(1:end-1)~=stimtypeorg(2:end); true]);
frame(2:2:end) = frameorg(novelend*2);
stimtype = stimtypeorg(novel);
%%%%%%% TO ANNA: commend following lines and uncommend following lines to get trackball info %%%
% load('trackball.mat')
% step=info.max_idx/(numel(tsamp)-1);
% velocity = interp1(0:step:info.max_idx,vsmooth,0:info.max_idx); 
%%%%%%%%%%%%%%%%
%% 
lastsignal=size(stim,2) -size(sig,1);
if lastsignal
    sig(size(sig,1):size(sig,1)+lastsignal,:)=0;
end

SNR = skewness(sig);
[~,I]= sort((SNR),'descend');
try
    sigplt_onepage(sig(:,I(1:30)),velocity,frame,stimtype); %
end
%% set paramenter for stimulus
prestim=ceil(CAframeHz*1); %use 1s baseline
stimON=median(frame(3:2:end)-frame(1:2:(end-2))); % StimON duration
seg=prestim+median(frame(5:4:end)-frame(1:4:(end-4))); % segment=prestim+info.frame(TTL ON+TTL off )

Var=numel(unique(stimtype));% calculate types of stimulus, orientations, contrast, etc
%rep=floor(min(2*numel(stimtype),numel(info.frame))/2/Var); %calculate repetitions
figure;
m=histogram(stimtype,.5:Var+.5);
title('Stim type')
rep = min(m.Values);
% rep=12;
sigT=zeros(seg*rep,Var,ncell); % build sigT based on signals from sig and organized according to segment
frame_on=zeros(Var,rep); %tag frame_on of each seg
win_sig=[prestim prestim+stimON-1];
%% reorganize stimulus according to stimulus type
last=frame(end-1)+seg-prestim;
if last>size(sig,1)
    disp(sprintf('padding zeros to the last %d frames', last-size(sig,1)));
    sig(end+1:last,:)=0;
    velocity(end+1:last)=0;

end

for j=1:Var
    ori=find(stimtype==j);
    on=ori*2-1;
    frame_on(j,:)=frame(on(1:rep)); % the frame marking the StimON timepoint size 1*rep
    temp=ones(seg,1)*frame_on(j,:)+ (0:seg-1)'*ones(1,rep)-prestim; % array seg*rep6
    sigT(:,j,:)=sig(temp,:);   %(seg*rep)*Var*ncell
end

sigT=reshape(sigT,seg,rep,Var,ncell);
%% syn with ball/eye motion
if isempty(velocity)
    run.matrix=[];
    run.thr = 0;
    run.data = [];
else
%     load([fn '_ball.mat'],'run');
%     if ~exist('run','var')
        run.data = [];
        for k=1:seg
            run.data = cat(3,run.data,velocity(frame_on'-prestim+k-1));%frame_on=zeros(Var,rep)
        end
        save([fn '_ball.mat'],'run','-append');
%     end
    [run.matrix,hrun,run.thr]=runplt2(run.data,win_sig);
end
%% parse sigF according to different types of info.var

if info.steps(2)>1
    [ sigF,variant,hparse]=parseF(sigT,win_sig,matrix);
else
    sigF=sigT; % sigF:seg,rep,Var,ncell
    variant=strcmp(info.var(1),'Orientation');
end

%% clean data
% sigF = cleandata(sigF,win_sig);

%% save data
folder=fullfile(fname(1:pos(end)-1),['picked' num2str(ncell)]);
if exist(folder,'dir')
    i = numel(dir ([folder '*']))+1;
    folder = [folder '_' num2str(i)];
end
mkdir(folder);

save(fullfile(folder,'peakSI.mat'),'sigF','run',...
    'prestim','stimON','win_sig','variant');

if exist('hsig','var')
    saveas(hsig,fullfile(folder,hsig.Name),'fig');
end

if exist('hrun','var')
    for i=1:numel(hrun)
        saveas(hrun(i),fullfile(folder,hrun(i).Name),'fig');
    end
end
if exist('hpick','var')
    saveas(hpick,fullfile(folder,hpick.Name),'fig');
end
if exist('newf','var')
    saveas(img, fullfile(folder,[newf '.png']));
    saveas(img, fullfile(folder,[newf '.fig']));
end
%% plotting
caanalysisplot(folder,find(gd1.f1));