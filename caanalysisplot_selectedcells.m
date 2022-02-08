 function caanalysisplot(varargin)
 % caanalysisplot(foldername, Cor, figskip)
 % input: foldername
 % input: Cor, selected good cells
 % input: figskip, not plotting each cell TC
 % example caanalysisplot('.',f1.pair0&f1.pair1&pair2,0)
if nargin<1
    [~,folder1]=uigetfile('peakSI.mat');
    [~,folder2]=uigetfile('peakSI.mat');
    rewrite=1;
else
    folder1 = varargin{1};
    folder2 = varargin{1};
    rewrite =0;
end
load(fullfile(folder1,'peakSI.mat'));
load(fullfile(folder2,'peakSI.mat'));

if nargin<2 
    Cor =1:size(sigF,4);% keep the last one even most times it's noise background
 % sigF: seg,rep,Var,ncell
else 
    Cor = varargin{2};
    if isempty(Cor) %
    Cor =1:size(sigF,4);
    end    
end  

if nargin<3
    figskip = 1;  % default to print all tuning 
else
    figskip = varargin{3};
end

if nargin<4
    plotmore = 0;
end
% sigF = sigF(:,:,:,Cor);
%% 
%save([fname '_' num2str(numel(Cor)) 'signals.mat'],'sigF','matrix','win_sig','Cor');
if ~exist('matrix','var') 
    matrix = run.matrix;
end

if ~exist('sigF','var')
    variant=1;
    try
        sigF = SI.sigF;
        SI = rmfield(SI,'sigF');
    catch
        [file,path]=uigetfile('*signals.mat','select sigF file');
        signals=load(fullfile(path,file))
        matrix = signals.matrix;
        win_sig = signals.window;
        sigF= signals.sigF;
    end
end
 
sigF(isnan(sigF))=0;

if win_sig(end) == size(sigF,1)
    temp=inputdlg(sprintf('segment=%d %d Manuel:',win_sig));
    win_sig = str2num(temp{1});
end
%% plot figures and calculate Orientation selection or Contrast invariant

wrunning = (~isempty(matrix)) & sum(matrix(:));
drawingoption=double(wrunning)*10+variant;
% variant=11;
% drawingoption=11;
% drawingoption=(~isempty(matrix)& sum(matrix(:)))*10+variant;

switch drawingoption
    case 0%'no running with contrast map'
        hsigF=sigFplt(sigF,[],win_sig,Cor,folder1);  % sigF:seg,rep,Var,ncell
        hsigF=sigFplt(sigF,[],win_sig,Cor,folder2); 
        [SI.baseline,SI.baselineSD,SI.peak,SI.error,cal_ER_option]=cal_ER(sigF,win_sig);%peak: 1*Var*ncell
        
        hpk(1)=fitplt(SI.peak,SI.error,Cor);
        hpk(1).Name='Contrast tuning under best orientation';
        
    case 1 %'no running with orientation map'
if ~figskip
         hsigF=sigFplt(sigF,[],win_sig,Cor,folder1);
         hsigF=sigFplt(sigF,[],win_sig,Cor,folder2);  % sigF:seg,rep,Var,ncell
end
       [SI.baseline,SI.baselineSD,SI.peak,SI.error,cal_ER_option]=cal_ER(sigF,win_sig);  %peak:1*Var*ncell
        [SI.OSI,SI.gOSI,SI.DSI,SI.gDSI,SI.pref_ori,SI.pref_dir]=calOSI(SI.peak);
                
        hpk(1)=polarplt(SI.peak,Cor);
        hpk(1).Name='Orientation tuning (polar)';
         
        hpk(2)=fitplt(SI.peak,SI.error,Cor);
        hpk(2).Name='Orientation tuning';

if plotmore
        hpk(3)=figure('Name','Preferred Orientation');
        h1 = histogram(SI.pref_ori);
        h1.FaceColor=[ 0.5 0.5 0.5];
        legend('anethtized','Location','north','Orientation','horizontal');
        xlabel('Orientation');ylabel('cell count');
        
        hpk(4)=figure('Name','Preferred Direction');
        h1 = histogram(SI.pref_dir);
        h1.FaceColor=[ 0.5 0.5 0.5];
        legend('anethtized','Location','north','Orientation','horizontal');
        xlabel('Direction');ylabel('cell count');
end
    case 2 %'no running with both contrast&orientation map'
if ~figskip
        hsigF=sigVarplt(sigF,info.steps,win_sig,matrix,Cor,folder1); 
        hsigF=sigVarplt(sigF,info.steps,win_sig,matrix,Cor,folder2);  % sigF:seg,rep,Var,ncell
end
        [SI.baseline,SI.baselineSD,SI.peak,SI.error,cal_ER_option]=cal_ER(sigF,win_sig);%peak:1*Var*ncell
%         
        hpk(1)=fitplt(SI.peak,SI.error,Cor);           %peak: 1*Var*ncell
        hpk(1).Name='Orientation tuning under diff contrast-running';
        
        
        [SI.OSI,SI.gOSI, SI.DSI,SI.gDSI,SI.pref_ori,SI.pref_dir]=calOSI(SI.peak);
        hpk(2)=polarplt(SI.peak,Cor);
        hpk(2).Name = 'Still only';
        hpk(3)=fitplt(SI.gOSI,[],Cor);
        hpk(3).Name= 'gOSI under diff contrast';
        
        hpk(4)=figure('Name','Preferred Orientation');
        h1 = histogram(SI.pref_ori);
        h1.FaceColor=[ 0.5 0.5 0.5];
        legend('anethtized','Location','north','Orientation','horizontal');
        xlabel('Orientation');ylabel('cell count');
    case 10  %'running with contrast'
        [SI.baselineR,SI.baselineSDR,SI.peakR,SI.errorR,...
            SI.baselineS,SI.baselineSDS,SI.peakS,SI.errorS]= sigFcmp(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
if ~figskip
        hsigF=sigFplt(sigF,matrix,win_sig,Cor,folder1);
        hsigF=sigFplt(sigF,matrix,win_sig,Cor,folder2);  % sigF:seg,rep,Var,ncell
end
%         SI.peakR=fixpeak(SI.peakR);
%         SI.peakS=fixpeak(SI.peakS);
        hpk(1)=figure('Name','Contrast tuning under best orientation-running');hold on;
        hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
        
    case 11  %'running with orientation'
        
        [SI.baselineR,SI.baselineSDR,SI.peakR,SI.errorR,...
            SI.baselineS,SI.baselineSDS,SI.peakS,SI.errorS,cal_ER_option]= sigFcmp(sigF,win_sig,matrix);  % sigR:seg,1,Var,ncell
       
        %         SI.peakR=fixpeak(SI.peakR);
        %         SI.peakS=fixpeak(SI.peakS);
%           %for Lab computer
        Day1=load('D:\Box\Enhancement-Andrea\1847\032917\000\data\1847_329_000_1\picked171_sp+ca\peakSI.mat');
        Day2=load('D:\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_1\picked455_sp+ca_3\peakSI.mat');
        
        p1=load('D:\Box\Enhancement-Andrea\1847\329&404\1847_329_000_1&1847_329_003_1&&1847_404_000_1&1847_404_001_1_plane1_selected142.mat');


%        
%         %for home computer
%         Day1=load('C:\Users\UCL053324\Box\Enhancement-Andrea\1847\032917\000\data\1847_329_000_1\picked171_sp+ca\peakSI.mat');
%         Day2=load('C:\Users\UCL053324\Box\Enhancement-Andrea\1847\040417\000\1847_404_000_1\picked455_sp+ca_3\peakSI.mat');
% 
%         p1=load('C:\Users\UCL053324\Box\Enhancement-Andrea\1847\329&404\1847_329_000_1&1847_329_003_1&&1847_404_000_1&1847_404_001_1_plane1_selected142.mat');
        compare2days(Day1,Day2);
    case 12 %'running with both orientation and contrast/eye/etc'
if ~figskip
        hsigF=sigVarplt(sigF,info.steps,win_sig,matrix,Cor,folder1);
        hsigF=sigVarplt(sigF,info.steps,win_sig,matrix,Cor,folder2);  % sigF:seg,rep,Var,ncell
end          
        [SI.baselineR,SI.baselineSdR,SI.peakR,SI.errorR,~,...
            SI.baselineS,SI.baselineSdS,SI.peakS, SI.errorS,~, cal_ER_option]= sigFcmp(sigF,win_sig,matrix);  % calculate peak from prestim:prestim+stimON-1
%         SI.peakR=reshape(SI.peakR,info.steps(2),info.steps(1),[]);
%         SI.errorR=reshape(SI.errorR,info.steps(2),info.steps(1),[]);
%         
%         SI.peakS=reshape(SI.peakS,info.steps(2),info.steps(1),[]);
%         SI.errorS=reshape(SI.errorS,info.steps(2),info.steps(1),[]);
%         
%         hpk(1)=fitplt(cat(1,SI.peakR,SI.peakS),cat(1,SI.errorR,SI.errorS),Cor);
%         hpk(1).Name='Orientation tuning under various condition-running';
%        
        [SI.OSI_R,SI.gOSI_R,SI.DSI_R,SI.gDSI_R,SI.pref_ori_R,SI.pref_dir_R]=calOSI(SI.peakR);
        %hpk(2).Name = 'Running condition';
        [SI.OSI_S,SI.gOSI_S,SI.DSI_S,SI.gDSI_S,SI.pref_ori_S,SI.pref_dir_S]=calOSI(SI.peakS);
        %hpk(3).Name = 'Still condition';

%         hpk(4)=fitplt(cat(1,SI.gOSI_R,SI.gOSI_S),[],Cor);
%         hpk(4).Name='gOSI under diff contrast and running state';
end
%% save/rewrite some figure & parameters
folder1 = fullfile(folder1,cal_ER_option);
folder2 = fullfile(folder2,cal_ER_option);
if ~exist(folder1)
    mkdir(folder1)
   ~exist(folder2)
    mkdir(folder2)
end

% if rewrite
%     baselinewin=win_sig(1)-5:win_sig(1);
%     baseline=nanmean(sigF(baselinewin,:,:,:));
%     seg=size(sigF,1);
%     sigF=sigF-repmat(baseline,[seg 1 1 1 ]);
%     threshold=2*nanstd(sigF(baselinewin,:,:,:));% threshold size: 1,1,Var,ncell
%     sigF(repmat(mean(sigF(win_sig(1):win_sig(end),:,:,:))<threshold,[seg 1 1 1 ]))=NaN; % set
%     [b,a] = butter(12,.2,'low');           % IIR filter design
%     sigF = filtfilt(b,a,sigF);
%     save(fullfile(folder,'peakSI.mat'),'SI','sigF','variant','Cor','matrix','win_sig');
    %     picked= dir('*pick2ed*.fig');
    %     pickedimg=picked.name(1:end-4);
    %     cd ../..;
    %     img=openfig([strtok(pickedimg,'pi') '.fig']);
    %     axesObjs = get(img, 'Children');  %axes handles
    %     dataObjs = get(axesObjs, 'Children');
    %     ncell =( numel( dataObjs)-1)/2;
    %     bad = setdiff(1:ncell,Cor);
    %     temp=[dataObjs(ncell+1-bad),dataObjs(2*ncell+1-bad)];
    %     delete(temp);
    %     title(strrep(pickedimg,'_','\_'),'FontSize',16);
    %     saveas(img, fullfile(folder,[pickedimg '.fig']));
    %     saveas(img, fullfile(folder,[pickedimg '.png']));
    %     caanalysisplot(folder);
% else
%%
% if ~exist(fullfile(folder,'peakSI.mat'))
%         save(fullfile(folder,'peakSI.mat'),'SI','sigF','variant','Cor','matrix','win_sig','-v7.3');
%     else
%         save(fullfile(folder,'peakSI.mat'),'SI','-append'); %need to update SI calculation JSUN
%     end
%     if exist('hpk','var')
%         for i=1:numel(hpk)
%             try
%         %         savefig(hpk(i),fullfile(folder,hpk(i).Name));
%                 saveas(hpk(i),fullfile(folder,hpk(i).Name),'png');
%                 end
%             end
%     end
%%  end
%close all
