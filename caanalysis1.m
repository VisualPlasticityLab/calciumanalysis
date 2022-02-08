function folder=caanalysis1(fname,fast)

% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal, movement, and external stimulus
% 2. interpolate ca signal, derive smoothed or spiking version
% 3. syn with movement (ball/ eye motion) and external stimulus internal state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% default not calculating spiking activity
if nargin<2 
    fast = 1;
end
%% 1. load calcium signal, movement, and external stimulus
if nargin ==0
    [fname,p]=uigetfile('*.signals','load calcium signal data .signals');
else
    [p,fnam,ext]=fileparts(fname);
    fname=[fnam ext];
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
velocity=preprocess(fn); %velocity has ball info and info has stimulus info
%% 2. interpolate ca signal, derive smoothed or spiking version

sig = sig1;
ncell=size(sig,2);
Cor = 1:ncell;
sig=sig(:, Cor);
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
    frame = floor((info.frame-multiplane)/nplane)+1;
    velocity = interp1(xx,velocity,x,'previous');
    nframes = floor((Nframes-multiplane)/nplane)+1;
end

switch fast
    case 0 %using deconv and MLspike to extract Calcium signal
        if ~exist('sigML','var') %using deconvolve to extract Calcium signal
            fprintf('MLspiking,size %d %d,expecting %d sec\n',size(sig),round(prod(size(sig))/600));tic
            sigML = nan(size(sig));
            sigs = nan(size(sig));
            sigsp = nan(size(sig));
            for ii=1:ncell
                try
                    sigML(:,ii)=MLspike(sig(:,ii));
                catch
                         fprintf('cell %d did not process',ii);
                end               
                if mod(ii,50)==1
                    fprintf('MLspiking %d out of %d cells in %d secs\n',ii,ncell,round(toc))
                end
            end
            fprintf('Done deconvolveCa and MLspiking, took %d secs\n',round(toc))
            save(fname,'sigML','-append');
        end
    case .5
        if ~exist('sigsp','var') %using deconvolve to extract Calcium signal
            fprintf('deconvlve...,size %d %d\n',size(sig));tic
            sigML = nan(size(sig));
            sigs = nan(size(sig));
            sigsp = nan(size(sig));
            for ii=1:ncell
                try
                    [sigs(:,ii),sigsp(:,ii)]=deconvolveCa(sig(:,ii));
                catch
                    fprintf('cell %d did not process',ii);
                end
                if mod(ii,50)==10
                    fprintf('deConvolveCa %d out of %d cells in %d secs\n',ii,ncell,round(toc))
                end
            end
            fprintf('Done deConvolveCa, took %d secs\n',round(toc))
            save(fname,'sigsp','sigs','-append');
        end
    case 1
        fprintf('skip spiking extraction\n')
        sigML = nan(size(sig));
        sigs = nan(size(sig));
        sigsp = nan(size(sig));
end
%% 3. syn with movement (ball/ eye motion) and external stimulus internal state
CAframeHz =info.resfreq/info.recordsPerBuffer/nplane;
prestim = floor(CAframeHz*1);
%Correct stimtype obtained with closeloop
stimtype = info.stimtype;
stimtypeorg = stimtype;
frameorg=frame;
grey = max(stimtype);
blank = 0;
if sum(stimtype == grey)>= numel(stimtype)/2-1 % this means half of the stimulus is gray stim, which is actually blank
        closeloop=1;
        stimtype([0;diff(stimtype)]~=0 & stimtype ==grey)=blank;     % if a stimtype is grey but the one before is different, then it's actually a blank
    while sum(stimtype(diff(stimtype)==0) ==grey) >=1 %two greys in a row
        % a grey stimtype is considered blank if: the one before is also grey, but the one 2 frames before is blank
        stimtype([1;1;stimtype(1:end-2)]==0& [1;stimtype(1:end-1)==grey] & stimtype ==grey)=blank;
    end
    firststim = find(stimtype==0,1);
    if firststim>2
        stimtype(1:firststim-2) = [];
        frame(1:2*(firststim-2)) =[];
    end
    stimON_each = frame(2:2:end)-frame(1:2:(end-1)); % StimON AND stimOFF duration 
%     run_each = velocity(
    seg_each = frame(4:4:end)-frame(1:4:(end-3)); % segment=prestim+info.frame(TTL ON+TTL off )
else
    %stimON=prctile(diff(info.frame),80); % StimON duration
     %StimON duration based on regular type
     closeloop=0;
    stimON_each = frame(2:2:end)-frame(1:2:(end-1));
    seg_each = frame(3:2:end)-frame(1:2:(end-2)); % segment=prestim+info.frame(TTL ON+TTL off )
end 
hh=figure('Position',[200 200 600 800],'Name','Stim');
subplot(3,1,1);hold on;
plot(frameorg(1:numel(stimtypeorg)*2),reshape(repmat(stimtypeorg',2,1),[],1),'o-');
plot(frame(1:numel(stimtype)*2),reshape(repmat(stimtype',2,1),[],1),'x-');
title(sprintf('Stimtype correction for closeloop experiment:%s',strrep(fn,'_','-')))
xlim([1 1000])
subplot(3,1,2);hold on;
histogram(stimON_each(1:2:end));histogram(seg_each);
legend('stimONduration','eachTrialduration')
stimON = floor(prctile(stimON_each(1:2:end),10));
seg = prestim + floor(prctile(seg_each,10));
Var=numel(unique(stimtype(stimtype>0)));% excluding blankstim(o) and calculate types of stimulus, orientations, contrast, etc
%rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions
subplot(3,1,3);hold on;
m=histogram(stimtype(stimON_each>stimON-CAframeHz/2),1:Var+1);
rep = min(m.Values);
sigF=zeros(seg*rep,Var,ncell); % build sigF based on signals from sig and organized according to segment
sigspF=zeros(seg*rep,Var,ncell); % build sigF based on signals from sig and organized according to segment
sigsF=zeros(seg*rep,Var,ncell); % build sigF based on signals from sig and organized according to segment
frame_on=zeros(Var,rep); %tag frame_on of each seg
win_sig=[prestim+1 prestim+stimON];

%% Plot example traces coupled with motion and stimulus
stim=zeros(1,nframes);
%for each stimtype, it sends signal to info.frame for its beginning and ending
for ii=1:(numel(frame)/2)
    stim(frame(ii*2-1):frame(ii*2))=stimtype(ii); 
end
if size(stim,2)> nframes
    stim(nframes+1:end)=[];
end
example=unique(randi(numel(Cor),1,50));

try
    sigplt_onepage(sig(:,example),velocity,frame,stimtype,Cor(:,example)); %
end
%% reorganize stimulus/ball motion according to stimulus type
last=frame(end-1)+seg-prestim;
if last>nframes
    disp(sprintf('padding zeros to the last %d frames', last-size(sigs,1)));
    sig(end+1:last,:)=0;
    sigs(end+1:last,:)=0;
    sigsp(end+1:last,:)=0;
    velocity(:,end+1:last)=0;
    sigML(end+1:last,:)=0;
end
% if  mod(numel(stimtype),2)
%     stimtype(end) = [];
% end
for j=1:Var
    ori=find(stimtype==j);
    % only keep the ones that has duration >=set ONduration
    if closeloop
        ori = ori(stimON_each(ori)>=stimON-CAframeHz/2);
%         tempvelocity = mean(velocity(info.frame(on)*ones(1,stimON)+ones(numel(on),1)*(0:stimON-1)),2); %rep*StimONduration
%         runtrial = tempvelocity>1;
%        frame_on(j,:)=info.frame(on(1:rep)); % the frame marking the stimON timepoint size 1*rep
    end
    on=ori*2-1;
    frame_on(j,:)=frame(on(1:rep)); % the frame marking the stimON timepoint size 1*rep
    temp=ones(seg,1)*frame_on(j,:)-prestim + (0:seg-1)'*ones(1,rep); % array seg*rep
    if temp(1,1)==0
        temp(1,1)=1;
    end
    sigO(:,j,:)=sig(temp,:);   %(seg*rep)*Var*ncell
    sigsO(:,j,:)=sigs(temp,:);
    sigspO(:,j,:)=sigsp(temp,:);
    sigML0(:,j,:)=sigML(temp,:);
end
sigO=reshape(sigO,seg,rep,Var,ncell);
sigsO=reshape(sigsO,seg,rep,Var,ncell);
sigspO=reshape(sigspO,seg,rep,Var,ncell);
sigMLO=reshape(sigML0,seg,rep,Var,ncell);

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
     [sigF,variant]=parseF(sigO(:,:,1:end-1,:),win_sig,run.matrix);
disp('define which version to show');
else
    sigF=sigO; % sigF:seg,rep,Var,ncell
    sigsF=sigsO; % sigF:seg,rep,Var,ncell
    sigspF=sigspO;
    sigMLF=sigMLO;
    variant=strcmp(info.var(1),'Orientation');
end

%% save data
folder=fullfile(fname(1:pos(end)-1),['picked' num2str(ncell)]);
if exist(folder,'dir')
    ii = numel(dir ([folder '*']))+1;
    folder = [folder '_' num2str(ii)];
end
mkdir(folder);

save(fullfile(folder,'peakSI.mat'),'sigF','run',...
    'prestim','stimON','win_sig','variant');

switch fast
    case .5
        save(fullfile(folder,'peakSI.mat'),'sigsF','sigspF','-append')
    case 0
        save(fullfile(folder,'peakSI.mat'),'sigMLF','-append')
end

saveas(hh,fullfile(folder,hh.Name),'fig');

if exist('hsig','var')
    saveas(hsig,fullfile(folder,hsig.Name),'fig');
end

if exist('hrun','var')
    for ii=1:numel(hrun)
        saveas(hrun(ii),fullfile(folder,hrun(ii).Name),'fig');
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
switch fast
    case 1
        caanalysisplot(folder);
    case .5
%          caanalysisplot(folder);
         caanalysisplot_sp(folder);
    case 0
%          caanalysisplot_sp(folder);
        caanalysisplot_ML(folder);
        % caanalysisplot_s(folder);
 end
