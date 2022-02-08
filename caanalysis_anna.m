function folder=caanalysis(fname)
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal
% 2. syn with external stimulus
% 3. syn with internal state: ball motion/ eye motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% load calcium signal and open recorded image
if nargin ==0
    [fname,p]=uigetfile('*.signals','load calcium signal data .signals');
else
    [p,fnam,ext]=fileparts(fname);
    fname=[fnam ext];
end
if ~isempty(p)
    cd(p)
end
load(fname,'-mat');
pos = strfind(strtok(fname,'c'),'_');
try
    fn = fname(1:pos(3)-1);
catch
    fn = strtok(fname,'.');
end
global info;
sbxread(fn,1,1);
movement=preprocess(fn);

global info;
if ~exist('sig','var')
    try
        sig=sig_chunk;
%         sigs=sigs_chunk;
%         sigsp=sigsp_chunk;
    catch
        sig = rtdata;
    end
end
%% cleaning data

if ~exist('sigs','var') 
sigsp = zeros(size(sig));
sigs = zeros(size(sig));
% sigML = zeros(size(sig));
sprintf('Smoothing and Deconvovling, expected %d secs',round(size(sig,2)*6))
tic
% for i=1:size(sig,2)
%     try
%     [sigs(:,i),sigsp(:,i)]=deconvolveCa(sig(:,i));
%     catch
%         sprintf('Cell%d did not processed',i);
%     end
% end
sprintf('Done Smoothing and Deconvovling, took %d secs',round(toc))
% for i=1:ncell
%     sigML(:,i)=MLspike(sig(:,i));
% end
% else
%     sigs=interp1(x,sigs,xx,'spline','extrap');
%     sigsp = interp1(x,sigsp,xx,'spline','extrap');
% end
% save(fname,'sigs','sigsp','-append');
end

%%
ncell=size(sigs,2);
Cor = 1:ncell;
sig=sig(:, Cor);

nframes = info.max_idx+1;
if info.volscan == 0 & nframes/size(sig,1) < 2
    x=1:size(sig,1);
else
%     nplane = info.otparam(3);
    multiplane = str2num(fname(pos(end)-1));
    nplane = round(nframes/size(sig,1));
    x=multiplane:nplane:nframes;
end
xx=1:nframes;
%     figure;hold on;plot(x,sig(:,67),'Linewidth',2,'Color',[ 0 0 1 .7])
%     plot(xx,sig(:,67),'Color',[ 0 0 0 .7])
sig=interp1(x,full(sig),xx,'spline','extrap');
sigs=interp1(x,sigs,xx,'spline','extrap');
sigsp = interp1(x,sigsp,xx,'spline','extrap');

%% sync movement with stimulus data and plot
CAframeHz =info.resfreq/info.recordsPerBuffer;
if movement 
    load([fn '_ball.mat'],'velocity');
    if ~exist('velocity','var')
        load([fn '_ball.mat'],'ball','time');
        speed=abs(ball(1:end-1))/192*2./diff(time)'; %size of the imaging area 192pixel, 2cm
        sp=conv(speed,ones(1,ceil(CAframeHz*3))/ceil(CAframeHz*3),'same'); %conv the speed into with a 1s filter
        if numel(sp) > nframes
            velocity=downsample(sp,2);
        else
            velocity=sp;
        end
    end
else
    velocity=[];
end

if size(velocity,2)> nframes
    velocity(nframes+1:end)=[];
end
%% sync stimulus with data  
prestim = floor(CAframeHz*1);

%Correct stimtype obtained with closeloop
info.stimtypeorg=info.stimtype;
info.frameorg=info.frame;

grey = max(info.stimtype);
blank = 0;
if sum(info.stimtype ==grey)>= numel(info.stimtype)/2-1 % this means blank is counted as grey
    % if a stimtype is grey but the one before is different, then it's actually
        % a blank 
        closeloop=1;
        info.stimtype([0;diff(info.stimtype)]~=0 &info.stimtype ==grey)=blank;
    while sum(info.stimtype(diff(info.stimtype)==0) ==grey) >=1 %two greys in a row
        % if a grey stimtype is considered blank if: the one before is also grey, but the one 2 frames before is blank
        info.stimtype([1;1;info.stimtype(1:end-2)]==0& [1;info.stimtype(1:end-1)==grey] & info.stimtype ==grey)=blank;
    end
    firststim = find(info.stimtype==0,1);
%     firstframe = info.frame(find(diff(info.frame)>CAframeHz,1)+1);
    if firststim>2
        info.stimtype(1:firststim-2) = [];
        info.frame(1:2*(firststim-2)) =[];
    end
 %StimON duration based on close-loop
    stimON_each = info.frame(2:2:end)-info.frame(1:2:(end-1)); % StimON AND stimOFF duration 
    run_each = velocity(
            tempvelocity = mean(velocity(info.frame(on)*ones(1,stimON)+ones(numel(on),1)*(0:stimON-1)),2); %rep*StimONduration

    seg_each = info.frame(4:4:end)-info.frame(1:4:(end-3)); % segment=prestim+info.frame(TTL ON+TTL off )
else
    %stimON=prctile(diff(info.frame),80); % StimON duration
     %StimON duration based on regular type
     closeloop=0;
    stimON_each = info.frame(2:2:end)-info.frame(1:2:(end-1));
    seg_each = info.frame(3:2:end)-info.frame(1:2:(end-2)); % segment=prestim+info.frame(TTL ON+TTL off )
end 
hh=figure('Position',[200 200 600 800],'Name','Stim');
subplot(3,1,1);hold on;
plot(info.frameorg(1:numel(info.stimtypeorg)*2),reshape(repmat(info.stimtypeorg',2,1),[],1),'o-');
plot(info.frame(1:numel(info.stimtype)*2),reshape(repmat(info.stimtype',2,1),[],1),'x-');
title(sprintf('Stimtype correction for closeloop experiment:%s',strrep(fn,'_','-')))
xlim([1 1000])
subplot(3,1,2);hold on;
histogram(stimON_each(1:2:end));histogram(seg_each);
legend('stimONduration','eachTrialduration')
stimON = floor(prctile(stimON_each(1:2:end),10));
seg = prestim + floor(prctile(seg_each,10));
Var=numel(unique(info.stimtype(info.stimtype>0)));% excluding blankstim(o) and calculate types of stimulus, orientations, contrast, etc
%rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions
subplot(3,1,3);hold on;
m=histogram(info.stimtype(stimON_each>stimON-CAframeHz/2),1:Var+1);
rep = min(m.Values);
sigF=zeros(seg*rep,Var,ncell); % build sigF based on signals from sig and organized according to segment
sigspF=zeros(seg*rep,Var,ncell); % build sigF based on signals from sig and organized according to segment
sigsF=zeros(seg*rep,Var,ncell); % build sigF based on signals from sig and organized according to segment
frame_on=zeros(Var,rep); %tag frame_on of each seg
win_sig=[prestim+1 prestim+stimON];

%% Plot example traces coupled with motion and stimulus
stim=zeros(1,nframes);
%for each stimtype, it sends signal to info.frame for its beginning and ending
for i=1:(numel(info.frame)/2)
    stim(info.frame(i*2-1):info.frame(i*2))=info.stimtype(i); 
end
if size(stim,2)> nframes
    stim(nframes+1:end)=[];
end
example=unique(randi(numel(Cor),1,10));
% hsig=figure('Name','Example');
sigplt(sig(:,example),[stim/max(stim);velocity/max(velocity)],Cor(example),sigs(:,example),sigsp(:,example)); %
%% reorganize stimulus/ball motion according to stimulus type
last=info.frame(end-1)+seg-prestim;
if last>nframes
    disp(sprintf('padding zeros to the last %d frames', last-size(sigs,1)));
    sig(end+1:last,:)=0;
    sigs(end+1:last,:)=0;
    sigsp(end+1:last,:)=0;
%     sigML(end+1:last,:)=0;
end
if  mod(numel(info.stimtype),2)
    info.stimtype(end) = [];
end
for j=1:Var
    ori=find(info.stimtype==j);
    % only keep the ones that has duration >=set ONduration
    if closeloop
        ori = ori(stimON_each(ori)>=stimON-CAframeHz/2);
        tempvelocity = mean(velocity(info.frame(on)*ones(1,stimON)+ones(numel(on),1)*(0:stimON-1)),2); %rep*StimONduration
        runtrial = tempvelocity>1;
        frame_on(j,:)=info.frame(on(1:rep)); % the frame marking the stimON timepoint size 1*rep
    else
        on=ori*2-1;
        frame_on(j,:)=info.frame(on(1:rep)); % the frame marking the stimON timepoint size 1*rep
    end
    temp=ones(seg,1)*frame_on(j,:)-prestim + (0:seg-1)'*ones(1,rep); % array seg*rep
    if temp(1,1)==0
        temp(1,1)=1;
    end
    sigO(:,j,:)=sig(temp,:);   %(seg*rep)*Var*ncell
    sigsO(:,j,:)=sigs(temp,:);
    sigspO(:,j,:)=sigsp(temp,:);
end
sigO=reshape(sigO,seg,rep,Var,ncell);
sigsO=reshape(sigsO,seg,rep,Var,ncell);
sigspO=reshape(sigspO,seg,rep,Var,ncell);

if ~movement
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
        save([fn '_ball.mat'],'run','velocity','-append');
%     end
    [run.matrix,hrun,run.thr]=runplt2(run.data,win_sig);
end
%% parse sigF according to different types of info.var
if info.steps(2)>1
%     [sigF,variant]=parseF(sigO(:,:,1:end-1,:),win_sig,run.matrix);
disp('define which version to show');
else
    sigF=sigO; % sigF:seg,rep,Var,ncell
    sigsF=sigsO; % sigF:seg,rep,Var,ncell
    sigspF=sigspO;
    variant=strcmp(info.var(1),'Orientation');
end
%% clean data
sigF = cleandata(sigF,win_sig);
%% save data
folder=fullfile(fname(1:pos(end)-1),['picked' num2str(ncell)]);
if exist(folder,'dir')
    i = numel(dir ([folder '*']))+1;
    folder = [folder '_' num2str(i)];
end
mkdir(folder);

save(fullfile(folder,'peakSI.mat'),'sigF','run',...
    'prestim','stimON','win_sig','variant');

% save(fullfile(folder,'peakSI.mat'),'sigF','sigsF','sigspF','run',...
%     'prestim','stimON','win_sig','variant');
saveas(hh,fullfile(folder,hh.Name),'fig');

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
caanalysisplot(folder);
% caanalysisplot_s(folder);
% caanalysisplot_sp(folder);
