% [fname,p] = uigetfile('*','images to make movie');
d=dir('green*movie.tif');
figure;
for k=1:numel(d)
    fname = d(k).name;
    %makeCalmovie(fname);
    make360movie(fname);
    saveas(gcf, [strtok(fname,'.') '_snapshot'], 'tiff')
end

