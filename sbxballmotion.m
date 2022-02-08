function ball = sbxballmotion(fn,sz)

load(fn,'-mat');            % should be a '*_ball.mat' file
if nargin>1
    data = data(sz{1},sz{2},:,:);
end
data = squeeze(data);
M = double(max(data(:)));
m = double(min(data(:)));
nframes = size(data,ndims(data));
ball = zeros(1,nframes);

% estimate threshold

th = zeros(1,100);
idx = round(linspace(1,size(data,3),100));
for(i=1:100)
    z = double(data(:,:,idx(i)));
    z = (z-m)/(M-m);
    th(i) = graythresh(z);
end
th = median(th);


% estimate translation...

z1 = double(data(:,:,1));
z1 = (z1-m)/(M-m);
z1 = uint8(im2bw(z1,th));

for(n=1:size(data,3)-1)
    z2 = double(data(:,:,n+1));
    z2 = (z2-m)/(M-m);
    z2 = uint8(im2bw(z2,th));
    [u,v] = fftalign(z1,z2);
    ball(n) = v+1i*u;
    z1 = z2;
end


% old method

% data = squeeze(data);
% ball = zeros(1,size(data,3));
% f_lap = fspecial('laplacian',.5); 
% 
% z1 = double(data(:,:,1));
% z1 = filter2(f_lap,z1,'valid');
% 
% for(n=1:size(data,3)-1)
%     n
%     z2 = double(data(:,:,n+1));
%     z2 = filter2(f_lap,z2,'valid');
%     [u,v] = fftalign(z1,z2);
%     ball(n) = v+1i*u;
%     z1 = z2;
% end
save(fn,'ball','-append');     % append the motion estimate data...
figure('Name',fn);plot(time,abs(ball))
display(['finished:' fn])