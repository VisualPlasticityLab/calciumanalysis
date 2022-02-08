%% make a small change
caanalysisplot('..',[],1);

%% read
d=dir('peakSI.mat');
if ~isempty(d)
    botheyepath=pwd;
else
    [~,botheyepath]=uigetfile('.mat','select PeakSI matfile for both eye');
end
variableInfo = who('-file',fullfile(botheyepath,'peakSI.mat'));
pos1 = strfind(botheyepath,filesep);

sigtype = 'sp';
load(fullfile(botheyepath,'peakSI.mat'),'matrix','win_sig','sigF','SI');
if isempty(matrix)
    matrix = logical(zeros(size(sigF,2),size(sigF,3))); %#ok
end

eye1_matrix = matrix(:,2:2:end);
eye2_matrix = matrix(:,1:2:end-1);

eye1_sigF = sigF(:,:,2:2:end,:); % from stimfile ipsieye, starting=2
eye2_sigF = sigF(:,:,1:2:end-1,:);% from stimfile contraeye, starting=1

%% calculate baseline using all running trials/ still trials
SIZ = size(eye1_sigF);
nframe = SIZ(1);
rep = SIZ(2);
nori = SIZ(3);
ncell = SIZ(4);

baselinewindow = 9:win_sig(1)-1; %baseline calculated from prestim 1s or2s depdening on code. right now its last 2sec of pre window (first second might have a post response from prior stim)
%  baselinewindow = win_sig(2)+1:nframe; %baseline calculated from poststim 2s
sigwindow = win_sig(1):win_sig(2);

eye1.trial = reshape(eye1_matrix,rep*nori,1);
eye2.trial = reshape(eye2_matrix,rep*nori,1);
eye1_sigBase = reshape(mean(eye1_sigF(baselinewindow,:,:,:)),rep*nori,ncell);
eye2_sigBase = reshape(mean(eye2_sigF(baselinewindow,:,:,:)),rep*nori,ncell);

eye1_BaseR = eye1_sigBase(eye1.trial==1,:);
eye1_BaseS = eye1_sigBase(eye1.trial~=1,:);
eye2_BaseR = eye2_sigBase(eye2.trial==1,:);
eye2_BaseS = eye2_sigBase(eye2.trial~=1,:);

eye1.base_R = mean(eye1_BaseR);
eye1.base_S = mean(eye1_BaseS);
eye1.baseSEM_R = std(eye1_BaseR)/sqrt(sum(eye1.trial==1));
eye1.baseSEM_S = std(eye1_BaseS)/sqrt(sum(eye1.trial==0));

eye2.base_R = mean(eye2_BaseR);
eye2.base_S = mean(eye2_BaseS);
eye2.baseSEM_R = std(eye2_BaseR)/sqrt(sum(eye2.trial==1));
eye2.baseSEM_S = std(eye2_BaseS)/sqrt(sum(eye2.trial==0));

%% get stimulus based on orientation as well as running/still status
eye1.peak_R = squeeze(SI.peakR(:,2:2:end,:));
eye2.peak_R = squeeze(SI.peakR(:,1:2:end-1,:));

eye1.peak_S = squeeze(SI.peakS(:,2:2:end,:));
eye2.peak_S = squeeze(SI.peakS(:,1:2:end-1,:));

eye1.peakSEM_R = squeeze(SI.SdR(:,2:2:end,:))./repmat(sqrt(sum(eye1_matrix==1))',1,ncell);
eye2.peakSEM_R = squeeze(SI.SdR(:,1:2:end-1,:))./repmat(sqrt(sum(eye2_matrix==1))',1,ncell);

eye1.peakSEM_S = squeeze(SI.SdS(:,2:2:end,:))./repmat(sqrt(sum(eye1_matrix==0))',1,ncell);
eye2.peakSEM_S = squeeze(SI.SdS(:,1:2:end-1,:))./repmat(sqrt(sum(eye1_matrix==0))',1,ncell);

%% excel output -1baseline mean with error, 1 stim mean with error (total 12 ), 
ipsi  = [ eye1.base_R; eye1.baseSEM_R; reshape([eye1.peak_R' ;eye1.peakSEM_R'],ncell,[])' ; eye1.base_S; eye1.baseSEM_S; reshape([eye1.peak_S' ;eye1.peakSEM_S'],ncell,[])' ];
contra = [ eye2.base_R; eye2.baseSEM_R; reshape([eye2.peak_R' ;eye2.peakSEM_R'],ncell,[])' ; eye2.base_S; eye2.baseSEM_S; reshape([eye2.peak_S';eye2.peakSEM_S'],ncell,[])' ];

% eye1
% eye2

cc_pre3 = cat(1,ipsi,contra);

%% compare using parametric

GOOD1 = (eye1.peak_R-eye1.peakSEM_R-eye1.base_R-eye1.baseSEM_R>0)|(eye1.peak_S-eye1.peakSEM_S-eye1.base_S-eye1.baseSEM_S>0);
GOOD2 = (eye2.peak_R-eye2.peakSEM_R-eye2.base_R-eye2.baseSEM_R>0)|(eye2.peak_S-eye2.peakSEM_S-eye2.base_S-eye2.baseSEM_S>0);

GOOD_p = (sum(GOOD1)+sum(GOOD2))>0;
find(GOOD_p)

'#above have significantly different distributions, below is ranksum '
%% compare using ranksum
cell_r= nan(nori,ncell);
for cell=1:ncell
    for ori=1:nori
    [~,eye1R]=ranksum(eye1_BaseR(:,cell), mean(eye1_sigF(sigwindow,eye1_matrix(:,ori)==1,ori,cell)),'tail','left','alpha',.05/nori/4);
    [~,eye1S]=ranksum(eye1_BaseS(:,cell), mean(eye1_sigF(sigwindow,eye1_matrix(:,ori)==0,ori,cell)),'tail','left','alpha',.05/nori/4);
    [~,eye2R]=ranksum(eye2_BaseR(:,cell), mean(eye2_sigF(sigwindow,eye2_matrix(:,ori)==1,ori,cell)),'tail','left','alpha',.05/nori/4);
    [~,eye2S]=ranksum(eye2_BaseS(:,cell), mean(eye2_sigF(sigwindow,eye2_matrix(:,ori)==0,ori,cell)),'tail','left','alpha',.05/nori/4);
    cell_r(ori,cell) = eye1R|eye1S|eye2R|eye2S;
    end
end
GOOD_r = sum(cell_r)>0;
find(GOOD_r)

%% all cells
% cd ..
% load('043magensp_all.mat')

%% tuning (OSI, DSI) for running and still
selec_cells = find(GOOD_p == 1); %pick output depending on parametric (GOOD_p) or ranksum (GOOD_r)
tunning_selection= [selec_cells', SI.OSI_S(:, (selec_cells))', SI.OSI_R(:, (selec_cells))', ...
                           SI.DSI_S(:, (selec_cells))', SI.DSI_R(:, (selec_cells))'];   
  
%% gOSI: [cell, cont-still-OSI, ipsi-still-OSI, cont-run-OSI, ipsi-run-OSI, cont-still-DSI, ...]   
globalOSI_S = (squeeze(SI.gOSI_S));
globalOSI_R = (squeeze(SI.gOSI_R));
globalDSI_S = (squeeze(SI.gDSI_S))';
globalDSI_R = (squeeze(SI.gDSI_R))';

gtunning_selection= [selec_cells',globalOSI_S((selec_cells),:), globalOSI_R((selec_cells),:),...
                                globalDSI_S((selec_cells),:), globalDSI_R((selec_cells),:)];

%% Peak response in the same format as above 
pref_orientation_still = ORI1(1,selec_cells)'; %get ori of selected cells
pref_orientation_run = ORI1(2,selec_cells)';

peak_response_still = squeeze(SI.peakS); %has contra first and ipsi second
peak_response_still_ipsi= peak_response_still((2:2:end),:); %just ipsi
peak_response_still_contra= peak_response_still((1:2:end),:); %just contra

peak_response_run = squeeze (SI.peakR);  %same as above
peak_response_run_ipsi= peak_response_run((2:2:end),:);
peak_response_run_contra= peak_response_run((1:2:end),:);

peak_select_run = [];
peak_select_still = [];
for ii = 1 : length(selec_cells)
    peak_select_run= [peak_select_run; peak_response_run_ipsi(pref_orientation_run(ii), selec_cells(ii)), ...
                                       peak_response_run_contra(pref_orientation_run(ii), selec_cells(ii))]; %#ok
    peak_select_still= [peak_select_still; peak_response_still_ipsi(pref_orientation_still(ii), selec_cells(ii)), ...
                                       peak_response_still_contra(pref_orientation_still(ii), selec_cells(ii))]; %#ok
end

%% ODI
 odi = [selec_cells', ODI1(:, selec_cells)']; % still/run