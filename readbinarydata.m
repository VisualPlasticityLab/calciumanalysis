ops = ops1{3};
[Ly Lx] = size(ops.mimg1);    % size of binary in x,y
nt0 = 0;
nimgbatch = 5000; % how many frames you can load at a time (will depend on RAM)
Nsteps = 20;
%%
figure;
clf;
fid  = fopen(ops.RegFile, 'r'); % opens the registered binary
while 1
  data = fread(fid, Ly*Lx*nimgbatch, '*int16');
  if isempty(data)
    break;
  end
  data = reshape(data, Ly, Lx, []);
  data = data(ops.yrange, ops.xrange, :);
  NT   = size(data,3);
  for j = 1:Nsteps:NT
    imagesc(data(:,:,j),[0 10000]);%,[1000 3000]);
    colormap gray
    title(sprintf('frame %d', j+nt0));
    axis square;
    drawnow;
    
    pause(0.25);
  end
  nt0 = nt0 + NT;
end
fclose all;