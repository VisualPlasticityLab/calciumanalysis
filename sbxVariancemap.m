function V=sbxVariancemap(fn)
sbxread(fn,1,1);        % read one frame to read the header of the image sequence
global info;            % this contains the information about the structure of the image
m=info.aligned.m;
Y=zeros(info.sz);
V=zeros(info.sz);
fac=sqrt(info.max_idx);
tic
for i=0:info.max_idx
    z =sbxread(fn,i,1);
    z= squeeze(z(1,:,:));
    Y = circshift(z,info.aligned.T(i+1,:)); % align the image
    temp=((double(Y-m))/fac).^2;
    V=V+temp;
end
save([fn '.align'],'V','-append');