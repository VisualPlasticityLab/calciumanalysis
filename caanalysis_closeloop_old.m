
%function ca_analysis;
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all;global info;
% sbxeyemotiondir;
% sbxballmotiondir;
% Gd = [2 9 3 10 4 11 18 65 16 19 17 23 24 31 32 30 39 59 50 51 53 54 62 40 60 47 41 42 43];

%% load calcium signal and open recorded image
[fanme,p]=uigetfile('.signals','load calcium signal data .signals');
load(fullfile(p,fanme),'-mat');


if ~exist('sig')
    sig=sig_chunk;
end
if size(sig,1)<size(sig,2)
    sig=sig';
end
total=size(sig,2);

if ~exist('Cor')
    Cor=1:total;
else
    bad=setdiff(1:total,Cor);
end

sig=sig(:,Cor);
%% sync with stimulus data and plot
try
    pos = strfind(strtok(fname,'c'),'_');
    fn = fname(1:pos(3)-1);
catch
    [temp,path]=uigetfile('*.sbx','Pick the corresponding sbx file');
    fname=fullfile(path,strtok(temp,'.'));
end
movement=preprocess(fname);
%%
CAframeHz =info.resfreq/info.recordsPerBuffer;
newtime=info.time(info.counting==1);
(info.frame(end)-info.frame(1))/CAframeHz - (newtime(end)-newtime(1))
CAlineHz = info.config.lines;
%just for plotting purpose,assign stimtype to all the ON frames
stim=zeros(1,info.max_idx+1); 


% [n1, m1] = size(info.frame);
% [n2, m2] = size(info.stimtype);
% n2 = n2*2;
% k = n1-n2;
% if k<0
%     info.stimtype(n1+1:n2) = [];
% else 
%     info.stimframe(n2+1:n1) = [];
% end

for i=1:(numel(info.frame)/2) 
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i);
end
%%%%%%% TO ANNA: comment following lines and uncomment the next few %%%
if movement
    try
        load([fname '.ball'],'-mat');
    catch
        load([fname '_ball'],'-mat');
        speed=abs(ball(1:end-1))/192*2./diff(time)'; %size of the imaging area 192pixel, 2cm
    end
    % load([fname '.eye'],'-mat');
    sp=conv(speed,ones(1,ceil(CAframeHz*1))/ceil(CAframeHz*1),'same'); %conv the speed into with a 1s filter
    if numel(sp) >= info.max_idx*2
        velocity=interp1(0:1/2:info.max_idx, sp(1:info.max_idx*2+1),0:info.max_idx); % downsampling the speed into the frame num
    else
        velocity=sp;
    end
else
    velocity=[];
end
%%%%%%%%%%% TO ANNA: uncommend following lines to get trackball info %%%
% load('trackball.mat')
% step=info.max_idx/(numel(tsamp)-1);
% velocity = interp1(0:step:info.max_idx,vsmooth,0:info.max_idx); 
%%%%%%%%%%%%%%%%
%% 
lastsignal=size(stim,2) -size(sig,1);
if lastsignal
    sig(size(sig,1):size(sig,1)+lastsignal,:)=0;
end

example=unique(randi(Cor(end),1,10));
hsig=sigplt(sig(:,example),[stim/max(stim);velocity/max(velocity)],Cor(example));

%hsig=sigplt(sig,[stim/max(stim);velocity/max(velocity)],Cor); %
% baseline=prctile(sig,20,1);
% sig=sig./repmat(baseline,size(sig,1),1)-1;

%% set paramenter for stimulus
ncell=size(sig,2);
prestim=ceil(CAframeHz*1); %use 1s baseline
stimON=median(info.frame(3:2:end)-info.frame(1:2:(end-2))); % StimON duration
seg=prestim+median(info.frame(5:4:end)-info.frame(1:4:(end-4))); % segment=prestim+info.frame(TTL ON+TTL off )

Var=numel(unique(info.stimtype));% calculate types of stimulus, orientations, contrast, etc
%rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions
m=histogram(info.stimtype,1:Var);
rep = min(m.Values);
% rep=12;
sigT=zeros(seg*rep,Var,ncell); % build sigT based on signals from sig and organized according to segment
frame_on=zeros(Var,rep); %tag frame_on of each seg
win_sig=[prestim prestim+stimON-1];
%% reorganize stimulus according to stimulus type
last=info.frame(end-1)+seg-prestim;
if last>size(sig,1)
    sig(end+1:last,:)=0;
    disp(sprintf('padding zeros to the last %d frames', last-size(sig,1)));
end

for j=1:Var
    ori=find(info.stimtype==j);
    on=ori*2-1;
    frame_on(j,:)=info.frame(on(1:rep)); % the frame marking the prestim timepoint size 1*rep
    temp=ones(seg,1)*frame_on(j,:)+ (0:seg-1)'*ones(1,rep)-prestim; % array seg*rep6
    sigT(:,j,:)=sig(temp,:);   %(seg*rep)*Var*ncell
end

sigT=reshape(sigT,seg,rep,Var,ncell);
%% select cells based on their peak responses, let's not do this here but save it later
% 
% if ~isempty(strfind(fanme,'memmap')) & isempty(strfind(fanme,'pick2'))
%     hcell=openfig([fanme(1:findstr(fanme,'cell')+3) '.fig']);
%     if ~exist('bad','var')
%         bad=pick(fanme);
%     end
%     [bad2,hpick]=pick2(sigT,win_sig,Cor);
%     bad=unique([bad bad2]);
%     if ~isempty(bad)
%         Cor = setdiff(1:ncell,bad);
%         ncell= numel(Cor);
%         newf=[fanme(1:end-8) 'pick2ed' num2str(ncell)];
%         
%         axesObjs = get(hcell, 'Children');  %axes handles
%         dataObjs = get(axesObjs, 'Children');
%         temp=dataObjs(total+1-bad);
%         delete(temp);
%         savefig(hcell,[newf '.fig']);
%         
%         sig=sig(:,Cor);
%         sigT=sigT(:,:,:,Cor);
%         save([newf '.signals'],'Cor');
%     end
% end

%% syn with ball/eye motion

if ~movement
    run.matrix=[];
    run.thr = 0;
    run.data = [];
else
    load([fn '_ball.mat'],'run');
    if ~exist('run','var')
        run.data = [];
        for k=1:seg
            run.data = cat(3,run.data,velocity(frame_on'-stimON+k-1));%frame_on=zeros(Var,rep)
        end
        save([fn '_ball.mat'],'run','-append');
    end
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
sigF = cleandata(sigF,win_sig);

%% save  data
save([fname '_' num2str(numel(Cor)) 'signals.mat'],'sigF','matrix','thr','win_sig','Cor');

mark=findstr(fname,'_');
date=fname(mark(1)+1:mark(2)-1);
serial=fname(mark(2)+1:end);
folder=fullfile([ date '-' serial ],['picked' num2str(numel(Cor))],['drawingoption' num2str(drawingoption)]);
if ~exist(folder,'dir')
    mkdir(folder);
end

save(fullfile(folder,'peakSI.mat'),'SI','Gd');

for i=1:numel(hsigF)
    saveas(hsigF(i),fullfile(folder,hsigF(i).Name),'png');
end
for i=1:numel(hpk)
    savefig(hpk(i),fullfile(folder,hpk(i).Name));
    saveas(hpk(i),fullfile(folder,hpk(i).Name),'png');
end

if exist('hrun')
    for i=1:numel(hrun)
        saveas(hrun(i),fullfile(folder,hrun(i).Name),'png');
    end
end

if exist('hparse')
    for i=1:numel(hparse)
        saveas(hparse(i),fullfile(folder,hparse(i).Name),'png');
    end
end
if exist('hpick')
    saveas(hpick,fullfile(folder,hpick.Name),'png');
end
close all;

