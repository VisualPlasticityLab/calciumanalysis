function V=sbxVariancemap2plane(fn)

sbxread(fn,1,1);        % read one frame to read the header of the image sequence
global info;
if info.volscan==1
    nplanes=info.otparam(3);
else
    nplanes=1;
end
V=zeros([info.sz nplanes]);
fac=double(max(info.aligned.m(:)));
Y=zeros(info.sz);


tic
for i=0:info.max_idx
    z =sbxread(fn,i,1);
    z= squeeze(z(1,:,:));
    Y = circshift(z,info.aligned.T(i+1,:)); % align the image
    k=mod(i,nplanes)+1;
    temp=((double(Y-info.aligned.m(:,:,k)))/fac).^2; 
    V(:,:,k)=V(:,:,k)+temp;
end
toc
save([fn '.align'],'V','-append');