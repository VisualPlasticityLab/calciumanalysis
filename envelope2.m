function [S_up,S_up_loc,S_up_f] = envelope2(S_norm,value)
%S_norm: the signal
% value : 0-100 , 99%- max bound, 1% min bound

if ~exist('value')
    value = 99;
elseif value >100
    value = 99.9;
elseif value<0
    value = 0.1;
end
    
segments = floor(numel(S_norm)/500);
remains = mod(numel(S_norm),segments);
if remains>=1
    S_norm(end-mod(numel(S_norm),segments)+1:end)=[];
end
% k = factor(numel(S_norm))
S_seg = reshape(S_norm,[],segments);

for jj=1:size(S_seg,2)
    S_up(jj) = prctile(S_seg(:,jj),value);
    [~,loc] = min(abs(S_seg(:,jj)-S_up(jj)));
    S_up_loc(jj)= loc+size(S_seg,1)*(jj-1);
end

% S_up_f=interp1(S_up_loc',S_up','spline');
S_up_f=interp1(S_up_loc,S_up,1:numel(S_norm)+remains,'pchip');