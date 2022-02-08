function folder=caanalysis(fname,fast)
% load calcium signal and sync with the ball motion/ pupil dilation data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. load calcium signal, movement, and external stimulus
% 2. interpolate ca signal, derive smoothed or spiking version
% 3. syn with movement (ball/ eye motion) and external stimulus internal state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clc;
%% 1. load calcium signal, movement, and external stimulus
if nargin ==0
    [fname,p]=uigetfile('*.signals','load calcium signal data .signals');
else
    [p,fnam,ext]=fileparts(fname);
    fname=[fnam ext];
end
display(fname);load(fname,'-mat');
pos = strfind(strtok(fname,'z'),'_');
try
    fn = fname(1:pos(3)-1);
catch
    fn = strtok(fname,'.');
end
global info;
sbxread(fn,1,1);
velocity=preprocess1(fn); %velocity has ball info and info has stimulus info
% %% default to use sigsp
% if nargin<2 
%     fast = 0.5;
% end
% switch fast
%     case 0 % sig is the calcium trace
%     case .5 % sig is deconvolved trace
%          if ~exist('sigsp','var')
%             sig=extraSP(sig,fname,fast); %pnenskit
%         else
%              sig = double(sigsp);% suite2p 
%          end
%         fn = [fn '_sp'];
%     case 1 %sig is deconvolved using MLFextraction
%         sig=extraSP(sig,fname,fast); % pick a way
%         fn = [fn '_sp2'];
% end
sig_raw = double(sig);
sig = double(sigsp);
%% 2. interpolate ca signal, derive smoothed or spiking version
if ~exist('sig','var')
sig = sig_chunk;
end
ncell=size(sig,2);
Cor = 1:ncell;
sig=sig(:, Cor);
Nframes = info.max_idx+1;
if info.volscan == 0 & Nframes/size(sig,1) < 2
%     x=1:size(sig,1);
        nframes = Nframes;
        nplane = 1;
        frame = info.frame;
else
%     nplane = info.otparam(3);
    multiplane = str2num(fname(pos(end)-1));
    nplane = round(Nframes/size(sig,1));
    x=multiplane:nplane:Nframes;
    xx=1:Nframes;
    % sig=interp1(x,full(sig),xx,'spline','extrap');
    frame = floor((info.frame-multiplane)/nplane)+1;
    nframes = floor((Nframes-multiplane)/nplane)+1;
     if ~isempty(velocity)
        velocity = interp1(xx,velocity,x,'previous');
    end
end

%% 3. syn with movement (ball/ eye motion) and external stimulus internal state
CAframeHz =info.resfreq/info.recordsPerBuffer/nplane;

%Correct stimtype obtained with closeloop
stimtypeorg = info.stimtype;
stimtype = stimtypeorg;
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
    stimON_each = frame(2:4:end)-frame(1:4:(end-1)); % StimON AND stimOFF duration 
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
if closeloop
    title(sprintf('Stimtype correction for closeloop experiment:%s',strrep(fn,'_','-')))
else
    title(sprintf('Not closeloop exp, no stimtype correction'))
end
xlim([1 1000])

%stimtime calculation
subplot(3,1,2);hold on;
histogram(stimON_each);histogram(seg_each);
legend('stimONduration','eachTrialduration')

if closeloop
    prestim = floor(1*CAframeHz); % forcely define prestim = 1s before
elseif floor(prctile(stimON_each,50))+1>=floor(prctile(seg_each,50)) % this means stimON is not recorded but just the whole section
          fprintf('Delay:%.2f,Duration:%.1f,WaitInterval:%.1f(s)\n',info.stimtime)
          stimON=ceil(floor(prctile(seg_each,50))*(info.stimtime(2)-abs(info.stimtime(1)))/(info.stimtime(2)+info.stimtime(3)));
          prestim = floor(prctile(seg_each,50))-stimON;
%         temp=inputdlg(sprintf('stimON duration=%d Manuel:',stimON));
%         stimON = str2num(temp{1});    
else
    stimON=prctile(stimON_each,50);
    prestim = floor(prctile(seg_each,50))-stimON;
end

if info.stimtime(1)>0   % this means the extra gray period happens before the stimON
    prestim = prestim+floor(info.stimtime(1)*CAframeHz);
end

seg = prestim+stimON +prestim; % now define seg to be prestim +stimON+prestim
win_sig=[prestim+1 prestim+stimON];
fprintf('adjustement: stimON=[%d %d],eachTrialduration=%d\n',win_sig,seg)

% repetition calculation
Var=numel(unique(stimtype(stimtype>0)));% excluding blankstim(o) and calculate types of stimulus, orientations, contrast, etc
%rep=floor(min(2*numel(info.stimtype),numel(info.frame))/2/Var); %calculate repetitions
subplot(3,1,3);hold on;
m=histogram(stimtype(stimON_each>stimON-CAframeHz/2),1:Var+1);
rep = min(m.Values);
frame_on=zeros(Var,rep); %tag frame_on of each seg

%% Plot example traces coupled with motion and stimulus
stim=zeros(1,nframes);
%for each stimtype, it sends signal to info.frame for its beginning and ending
for ii=1:(numel(frame)/2)
    stim(frame(ii*2-1):frame(ii*2))=stimtype(ii); 
end
if size(stim,2)> nframes
    stim(nframes+1:end)=[];
end
example=unique(randi(numel(Cor)-2,1,10));

try
    sigplt_onepage(sig_raw(:,example),velocity,frame,stimtype,Cor(example)); %
saveas(gcf,[fn '_example_sig_raw.fig'])
saveas(gcf,[fn '_example_sig_raw.png'])

end
%% reorganize stimulus/ball motion according to stimulus type
last=frame(end-1)+seg-prestim;
npadding = last-size(sig,1);
% if last>nframes
if npadding >0
    fprintf('padding to the last %d frames\n', npadding );
    sig(end+1:last,:)=repmat(sig(end,:),npadding,1);
    sig_raw(end+1:last,:)=repmat(sig_raw(end,:),npadding,1);
%     velocity(:,end+1:last)=repmat(velocity(:,end),npadding,1);
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
    if temp(1,1)<1 % not enough time for the first prestim
        temp(temp(:,1)<1,1)=1;
    end
    sigO(:,j,:)=sig(temp,:);   %(seg*rep)*Var*ncell
    sigO_raw(:,j,:)=sig_raw(temp,:);   %(seg*rep)*Var*ncell

end
sig1=reshape(sigO,seg,rep,Var,ncell);
sig1_raw=reshape(sigO_raw,seg,rep,Var,ncell);

if isempty(velocity)
    run.matrix= [];
    run.thr = 0;
    run.data = [];
else
%     load([fn '_ball.mat'],'run');
%     if ~exist('run','var')  
        run.data = [];
        for k=1:seg
            timestamp = frame_on'-prestim+k-1;
            timestamp(timestamp<1)=1;
            run.data = cat(3,run.data,velocity(timestamp));%frame_on=zeros(Var,rep)
        end
%         save([fn '_ball.fmat'],'run','-append');
%     end
    [run.matrix,hrun,run.thr]=runplt2(run.data,win_sig);
end
%% parse sigF according to different types of info.var
if info.steps(2)>1
     [sigF,variant]=parseF(sig1,win_sig,run.matrix);
     sigF_raw=parseF(sig1_raw,win_sig,run.matrix);
else
    sigF=sig1; % sigF:seg,rep,Var,ncell
    sigF_raw=sig1_raw; % sigF:seg,rep,Var,ncell
    variant=strcmp(info.var(1),'Orientation');
end
if isempty(run.matrix)
    run.matrix = logical(zeros(size(sigF,2),size(sigF,3)));
end

%% save data
folder=fullfile(fname(1:pos(end)-1),['picked' num2str(ncell)]);
% switch fast
%     case 1
%     
%     case .5
%         folder = [folder '_sp'];
%     case 0
%         folder = [folder '_sp2'];
% end
    
        folder = [folder '_sp+ca'];

if exist(folder,'dir')
    ii = numel(dir ([folder '*']))+1;
    folder = [folder '_' num2str(ii)];
end
mkdir(folder);

save(fullfile(folder,'peakSI.mat'),'sigF','sigF_raw','run',...
    'prestim','stimON','win_sig','variant','info');
    

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
 %close all;
%% plotting

         caanalysisplot(folder,example,0); %close all;
         %caanalysisplot(folder,1:ncell,0); close all;

% switch fast
%     case 1
% %        [~,values]= AutoCellSelection(folder);
% %         caanalysisplot(folder,find(values(:,1)>.15&values(:,2)<.2)',0);
%         caanalysisplot(folder,example,0);
%     case .5
% %          caanalysisplot(folder);
%          caanalysisplot_sp(folder);
%     case 0
% %          caanalysisplot_sp(folder);
%         caanalysisplot_ML(folder);
%         % caanalysisplot_s(folder);
%  end