
% function folder=caanalysis(fname)
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load calcium signal and open recorded image
% if nargin ==0
[fname,p]=uigetfile('*.signals','load calcium signal data .signals');
% else
% [p,fname,ext]=fileparts(fname);
% fname=[fname ext];
% end
cd(p);
global info;

load(fname,'-mat');
pos = strfind(strtok(fname,'c'),'_');
fn = fname(1:pos(3)-1);
figf=[strtok(fname,'.') '.fig'];

if ~exist('sig')
    try
        sig=sig_chunk;
    catch
        sig = rtdata;
    end
end
if size(sig,1)<size(sig,2)
    sig=sig';
end
% sig=sig(:,Cor);
sig=sig(:,1:end-1);
ncell=size(sig,2);
nframes=size(sig,1);
%% sync with stimulus data and plot
movement=preprocess(fn);
CAframeHz =info.resfreq/info.recordsPerBuffer;

if movement
    load([fn '_ball.mat'],'ball','time','-mat');
    speed=abs(ball(1:end-1))/192*2./diff(time)'; %size of the imaging area 192pixel, 2cm
    sp=conv(speed,ones(1,ceil(CAframeHz*3))/ceil(CAframeHz*3),'same'); %conv the speed into with a 1s filter
    if numel(sp) >= info.max_idx*2
        velocity=downsample(sp,2);
    elseif multiplane
        velocity=interp1(1:(info.max_idx+1), sp,1:2:(info.max_idx+1));
    else
        velocity=sp;
    end
else
    velocity=[];
end
% 
% if ~exist('sigs')
%     sigsp=zeros(size(sig));
%     sigs=zeros(size(sig));
% for i=1:ncell
%     [sigs(:,i),sigsp(:,i)]=deconvolveCa(sig(:,i));
% end
% end
% for i=1:ncell
% sigs(:,i)=MLspike(sig(:,i));
% end
% sigsp=sigs;
sigs=sig;
if info.volscan == 0
    multiplane = 0;
else
    nplane = info.otparam(3);
    multiplane = str2num(fname(pos(end)-1));
    x=multiplane:nplane:info.max_idx+1;
    xx=1:info.max_idx+1;
    figure;hold on;plot(x,sigs(:,67),'Linewidth',2,'Color',[ 0 0 1 .7])
    sigs=interp1(x,sigs,xx,'spline');
    plot(xx,sigs(:,67),'Color',[ 1 0 0 .7])
end
if size(stim,2)> size(sigs,1)
    stim(size(sigs,1)+1:end)=[];
end
if size(velocity,2)> size(sigs,1)
    velocity(size(sigs,1)+1:end)=[];
end
if size(velocity,2)> size(stim,2)
    velocity=velocity(1:size(stim,2));
end
%% cleaning data
% baseline=prctile(sig,5,1);
% sig=sig./repmat(baseline,size(sig,1),1)-1;
% startfram = info.frame(1:2:end)';
% endfram = info.frame(2:2:end)';
% for i=1:ncell
%     sig(:,i)=func_preprocess_data(sig(:,i),startfram,endfram);   
% end
example=unique(randi(numel(Cor),1,10));
hsig=sigplt(sig(:,example),[stim/max(stim);velocity/max(velocity)],Cor(example),sigs(:,example),sigsp(:,example)); %

%% set paramenter for stimulus
prestim=ceil(CAframeHz*1); %use 1s baseline
stimON=info.frame(2)-info.frame(1); % StimON duration
seg = prestim+min(info.frame(3:2:end)-info.frame(1:2:(end-2))); % segment=prestim+info.frame(TTL ON+TTL off )
% stimON = ceil(CAframeHz*2);
% seg = prestim + stimON;
Var=numel(unique(info.stimtype));% calculate types of stimulus, orientations, contrast, etc
rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions
sigF=zeros(seg*rep,Var,ncell); % build sigF based on signals from sig and organized according to segment
% run = zeros(seg*rep,Var);
frame_on=zeros(Var,rep); %tag frame_on of each seg
win_sig=prestim:prestim+stimON;
%% reorganize stimulus according to stimulus type
last=info.frame(end-1)+seg-prestim;
if last>size(sigs,1)
    disp(sprintf('padding zeros to the last %d frames', last-size(sigs,1)));
    sigs(end+1:last,:)=0;
    sigsp(end+1:last,:)=0'
end

sigO=[];
for j=1:Var
    ori=find(info.stimtype==j);
    on=ori*2-1;
    frame_on(j,:)=info.frame(on(1:rep)); % the frame marking the prestim timepoint size 1*rep
    temp=ones(seg,1)*frame_on(j,:)+ (0:seg-1)'*ones(1,rep)-prestim; % array seg*rep6
    if temp(1,1)==0
        temp(1,1)=1;
    end
    sigO(:,j,:)=sigs(temp,:);   %(seg*rep)*Var*ncell
    sigspO(:,j,:)=sigsp(temp,:);
    %     run(j,:)= velocity(temp);
end
% run = reshape(run,seg,rep,Var);
sigO=reshape(sigO,seg,rep,Var,ncell);
sigspO=reshape(sigspO,seg,rep,Var,ncell);
%% syn with ball/eye motion
if ~movement
    run.matrix=[];
    run.thr = 0;
    run.data = [];
else
    variableInfo = who('-file', [fn '_ball.mat']);
%     if ismember('run', variableInfo)
%         load([fn '_ball.mat'],'run');
%     else
        run.adjusted=[];%frame_on=zeros(Var,rep)
        run.data = [];
        for k=1:seg
            run.data = cat(3,run.data,velocity(frame_on'+k-1));
        end
        %define an adjusted factor
        x=1:stimON;
        %y=exp(-x/20);
        y=exp(-(x-(stimON+1)/2).^2/(1+stimON)^2*2);
        factor=y./sum(y)*mean(x);
        for k=1:stimON
            run.adjusted=cat(3,run.adjusted,velocity(frame_on'+k-1)*factor(k));% run:rep,Var,stimON
        end
        [run.matrix,hrun,run.thr]=runplt2(run.data,run.adjusted);
        %[run.matrix,hrun]=runplt2(run);
         save([fn '_ball.mat'],'run','-append');
%     end
end
%% parse sigF according to different types of info.var
if info.steps(2)>1
    [ sigF,variant]=parseF(sigO(:,:,1:end-1,:),win_sig,run.matrix);
else
    sigF=sigO; % sigF:seg,rep,Var,ncell
    sigspF=sigspO;
    variant=strcmp(info.var(1),'Orientation');
end

baselinewin=win_sig(1)-5:win_sig(1)+5;
baseline=nanmin(sigF(baselinewin,:,:,:));
seg=size(sigF,1);
sigF=sigF-repmat(baseline,[seg 1 1 1 ]);
noisewin=win_sig(1)-5:win_sig(1)-1;
threshold=2*nanstd(sigF(noisewin,:,:,:));% threshold size: 1,1,Var,ncell
sigF(repmat(mean(sigF(win_sig(1):win_sig(end),:,:,:))<threshold,[seg 1 1 1 ]))=0; % set

%% save data

folder=fullfile(fname(1:pos(end)-1),['picked' num2str(ncell)]);

if exist(folder,'dir')
    folder = [folder '_2'];
end
mkdir(folder);

save(fullfile(folder,'peakSI.mat'),'sigF','sigO','sigspF','run',...
    'prestim','stimON','win_sig','variant');


if exist('hrun','var')
    for i=1:numel(hrun)
        saveas(hrun(i),fullfile(folder,hrun(i).Name),'png');
    end
end
if exist('hpick','var')
    saveas(hpick,fullfile(folder,hpick.Name),'png');
end
if exist('newf','var')
    saveas(img, fullfile(folder,[newf '.png']));
    saveas(img, fullfile(folder,[newf '.fig']));
end
%% plotting
caanalysisplot(folder);
