d=dir('peakSI.mat');
if ~isempty(d)
    botheyepath=pwd;
else
    [~,botheyepath]=uigetfile('.mat','select PeakSI matfile for both eye');
end
variableInfo = who('-file',fullfile(botheyepath,'peakSI.mat'));
pos1 = strfind(botheyepath,filesep);

sigtype = 'sp';
load(fullfile(botheyepath,'peakSI.mat'),'matrix','win_sig','sigF');
eye1.sigF = sigF(:,:,2:2:end,:); % from stimfile ipsieye, starting=2
eye2.sigF = sigF(:,:,1:2:end,:);% from stimfile contraeye, starting=1


if isempty(matrix)
    matrix = logical(zeros(size(sigF,2),size(sigF,3)));
end

eye1.matrix = matrix(:,2:2:end);
eye2.matrix = matrix(:,1:2:end-1);
%%

SIZ = size(eye1.sigF);
nframe = SIZ(1);
rep = SIZ(2);
nori = SIZ(3);
ncell = SIZ(4);


eye1.trial = reshape(eye1.matrix,rep*nori,1);
eye2.trial = reshape(eye2.matrix,rep*nori,1);
eye1.sigBase = reshape(mean(eye1.sigF(1:win_sig(1)-1,:,:,:)),rep*nori,ncell);
eye1.sigStim = reshape(mean(eye1.sigF(win_sig(1):win_sig(2),:,:,:)),rep*nori,ncell);
eye2.sigBase = reshape(mean(eye2.sigF(1:win_sig(1)-1,:,:,:)),rep*nori,ncell);
eye2.sigStim = reshape(mean(eye2.sigF(win_sig(1):win_sig(2),:,:,:)),rep*nori,ncell);

eye1.sigBase_R = eye1.sigBase(eye1.trial==1,:);
eye1.sigStim_R = eye1.sigStim(eye1.trial==1,:);
eye2.sigBase_R = eye2.sigBase(eye2.trial==1,:);
eye2.sigStim_R = eye2.sigStim(eye2.trial==1,:);

eye1.sigBase_S = eye1.sigBase(eye1.trial~=1,:);
eye1.sigStim_S = eye1.sigStim(eye1.trial~=1,:);
eye2.sigBase_S = eye2.sigBase(eye2.trial~=1,:);
eye2.sigStim_S = eye2.sigStim(eye2.trial~=1,:);

%aa = [eye1.sigBase_R
