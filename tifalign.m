function tifalign(fname,T,idx)

if ~exist('idx','var')
    info = imfinfo(fname);
    num_images = numel(info);
    idx = 1:num_images;
end
fname_new = sprintf('%s_%d-%d_%s',fname(1:end-4),idx(1),idx(end),'aligned');

if isempty(T)
    [m,T]= tifalignx(fname,idx);
    save([fname_new '.mat'],'m','T');
end

for i=idx
    k=imread(fname,i);
    img=circshift(k,T(idx==i,:));
    if i==idx(1)
        imwrite(img,[fname_new '.tif'],'tif');
    else
        imwrite(img,[fname_new '.tif'],'tif','writemode','append');
    end
end

