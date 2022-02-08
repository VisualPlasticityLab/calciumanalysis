function peak=fixpeak(peak)

%peak size: type(contrast/running)*ori*ncell
ncell = size(peak,3);
Var = size(peak,2);

[m_ctr,m_ori]=find(isnan(peak(:,:,1)));
for k=1:numel(m_ori)
    before=m_ori(k)-1;before=mod(before-1,Var)+1;
    after=m_ori(k)+1;after=mod(after-1,Var)+1;
    peak(m_ctr(k),m_ori(k),:)= (peak(m_ctr(k),before,:)+peak(m_ctr(k),after,:))/2;
end

%[pref0,pref_ori0]=nanmax(peak,[],2);
%peak2 =peak;
% figure
% ncolum=ceil(sqrt(ncell));
% for i=1:ncell
%     subplot(ncolum,ncolum,i);
%     peak2(:,1:end-1,i) = fit2peakgaussion(pref_ori0(1,:,i),squeeze(peak(1,1:end-1,i)));
% end