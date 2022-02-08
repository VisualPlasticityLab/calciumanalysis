function [eye1,eye2]=loadsigFdata

[~,ipsieyepath]=uigetfile('*.mat','select trace for ipsieye');


eye1=load(fullfile(ipsieyepath,'peakSI.mat'),'run','matrix','window','win_sig','sigF');
if isfield(eye1,'win_sig')
    win_sig1 = eye1.win_sig;
else
    load(fullfile(ipsieyepath,'peakSI.mat'),'window');
    win_sig1 = [eye1.window(1) eye1.window(end)];
end
if ~isfield(eye1,'matrix')
    load(fullfile(ipsieyepath,'peakSI.mat'),'run');
    eye1.matrix=run.matrix;
end
% try
%     [~,ipsieyepath2]=uigetfile('*.mat','Load more traces for ipsieye?');
%     eye1extra=load(fullfile(ipsieyepath2,'peakSI.mat'),'sigF','matrix','win_sig');
%     framestocapture=min(size(eye1extra.sigF,1),size(eye1.sigF,1));
%     eye1.sigF = cat(2, eye1.sigF(1:framestocapture,:,:,:),eye1extra.sigF(1:framestocapture,:,:,:));
%     if ~isfield(eye1extra,'matrix')
%         load(fullfile(ipsieyepath2,'peakSI.mat'),'run');
%         eye1extra.matrix=run.matrix;
%     end
%     eye1.matrix= cat(1,eye1.matrix,eye1extra.matrix);
% end

[~,contraeyepath]=uigetfile('*.mat','select trace for contraeye');
eye2=load(fullfile(contraeyepath,'peakSI.mat'),'matrix','window','win_sig','sigF');
if isfield(eye2,'win_sig')
    win_sig2 = eye2.win_sig;
else
    load(fullfile(contraeyepath,'peakSI.mat'),'window');
    win_sig2 = [eye2.window(1) eye2.window(end)];
end
if ~isfield(eye2,'matrix')
    load(fullfile(contraeyepath,'peakSI.mat'),'run');
    eye2.matrix=run.matrix;
end
% try
%     [~,contraeyepath2]=uigetfile('*.mat','Load more traces for contraeye?');
%     eye2extra=load(fullfile(contraeyepath2,'peakSI.mat'),'sigF','matrix');
%     framestocapture=min(size(eye2extra.sigF,1),size(eye2.sigF,1));
%     eye2.sigF = cat(2, eye2.sigF(1:framestocapture,:,:,:),eye2extra.sigF(1:framestocapture,:,:,:));
%     if ~isfield(eye2extra,'matrix')
%         load(fullfile(contraeyepath2,'peakSI.mat'),'run');
%         eye2extra.matrix=run.matrix;
%     end
%     eye2.matrix= cat(1,eye2.matrix,eye2extra.matrix);
% end

if isempty(eye1.matrix) | isempty(eye2.matrix)
    eye1.matrix = logical(ones(rep,nSteps));
    eye2.matrix = logical(ones(rep,nSteps));
end


[eye1.peakR,eye1.sigR,eye1.errorR,eye1.peakS,eye1.sigS,eye1.errorS]= sigFcmp(eye1.sigF,eye1.win_sig,eye1.matrix);
[eye2.peakR,eye2.sigR,eye2.errorR,eye2.peakS,eye2.sigS,eye2.errorS]= sigFcmp(eye2.sigF,eye2.win_sig,eye2.matrix);
[eye1.peak,eye1.sig,eye1.error]=cal_spER(eye1.sigF,eye1.win_sig);  %peak:1*Var*ncell
[eye2.peak,eye2.sig,eye2.error]=cal_spER(eye2.sigF,eye2.win_sig);  %peak:1*Var*ncell