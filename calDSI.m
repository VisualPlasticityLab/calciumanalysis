function DSI=calDSI(peak)
%peak size: type(contrast/running)*ori*ncell

dim=ndims(peak);
ncell=size(peak,dim);
ori=size(peak,dim-1);
type=size(peak,dim-2);

if ~(dim==3)
       peak=reshape(peak,[],ori,ncell);
end
    
%peak contrast*ori*ncell
if mod(ori,2)   %if with a blank stim
    ori=ori-1;
    peak = peak(:,1:ori,:);
end

[pref,pref_dir]=max(max(peak),[],2);
oppo_dir=mod(pref_dir+ori/2-1,ori)+1;

% A = pref;
% B = peak(:,sub2ind([ori,ncell],oppo_dir,1:ncell))
for i=1:ncell
    A=peak(:,pref_dir(i),i);
    B=peak(:,oppo_dir(i),i);
    DSI(:,i)=(A-B)./(A+B);
end



