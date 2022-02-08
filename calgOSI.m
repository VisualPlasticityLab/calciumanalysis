function gOSI=calOSI(peak)
%peak size: type(contrast/running)*ori*ncell
dim=ndims(peak);
ncell=size(peak,dim);
ori=size(peak,dim-1);
type=size(peak,dim-2);

if ~(dim==3)
       peak=reshape(peak,[],ori,ncell);
end

if mod(ori,2)   %if there is blank stimulus
    ori=ori-1;
    peak = peak(:,1:ori,:);
end

%peak contrast*ori*ncell

anglesRads = 2*pi*(1:ori)'/ori;
vector=exp(2i*anglesRads);
peak_p=permute(peak,[3 1 2]); %peak_p:ncell*contrast*ori
peak_p=reshape(peak_p,ncell*type,ori);
gOSI_p= abs(peak_p*vector)./sum(peak_p,2);
gOSI_p(gOSI_p==NaN)=0;
%gOSI_p= abs(peak_p*vector)./sum(peak_p,2).*max(peak_p,[],2);
gOSI_p=reshape(gOSI_p,1,ncell,type);
gOSI=permute(gOSI_p,[1 3 2]);





