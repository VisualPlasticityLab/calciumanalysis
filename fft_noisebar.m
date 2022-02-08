addpath(genpath('/home/idl/Documents/MATLAB'));

cropflag = 1; %0 for no cropping and 1 for cropping
GIFflag = 0; %1 for making animated figure
[FileName, PathName] = uigetfile('/data/*.mj2');

%% Cropping the video
[xgrid, ygrid] = crop_video(PathName, FileName, cropflag);

%% Loading the image
close all
vid = VideoReader([PathName, FileName]);
info = get(vid);
nFrames = info.Duration * info.FrameRate;
video = zeros(length(ygrid), length(xgrid), nFrames);

tic
count = 1;
while count <= nFrames
    vtemp  = readFrame(vid);
    video(:,:,count) = vtemp(ygrid, xgrid);
    if ~mod(count, 1000)
        disp(['Frame done: %', num2str(100*count/nFrames)])
    end
    count=count+1;
end
video = single(video);
disp(['Loading video and cropping took: ', num2str(toc)]);

%% set parameters (first try out rigid motion correction)
options_rigid = NoRMCorreSetParms('d1',size(video,1),'d2',size(video,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);

%% perform motion correction
[adjvid, shifts1, template1, options_rigid] = normcorre(video,options_rigid);
adjvid = single(adjvid);
clear video vtemp

%% Reorganizing the time and phase arrays
tic
load([PathName, FileName(1:end-4), '.mat'])
[~, ia, ~] = unique(frame_lst(:,1));
frame_lst = frame_lst(ia, :);
cinds = find(diff(frame_lst(:, 4)) < -180) + 1; %%finding the start of the cycles
bg_avg = mean(adjvid(:, :, cinds(end):end), 3); %%background (last 10 sec)
frame_lst = frame_lst(cinds(1):cinds(end)-1, :);
adjvid = adjvid(:, :, cinds(1):cinds(end)-1);
ext_frms = find(diff(frame_lst(:,1))==2);
adjvid(:, :, ext_frms) = [];
if size(adjvid, 3) < size(frame_lst, 1)
    disp('Frames are missing. Either repeat the experiment or CRY!')
end
disp(['Phase-time cleaning took: ', num2str(toc)]);

%% Doing F/Fbar
tic
FoverFbar = zeros(size(adjvid));
for ii = 1: size(adjvid, 3)
    FoverFbar(:, :, ii) = (adjvid(:, :, ii) - bg_avg)./bg_avg;
end
FoverFbar = single(FoverFbar);
clear avg_frames count ia 
disp(['Normalization took: ', num2str(toc)]);

%% Generating a GIF file
if GIFflag == 1
    tic
    making_GIF_file(PathName, FileName(1:end-4), FoverFbar)
    disp(['GIFation(!) took: ', num2str(toc)]);
end

%% Averaging frames over phase
tic
if size(frame_lst, 1) > size(FoverFbar, 3)
    frame_lst(end, :) = [];
    disp('Frame missing from frame_lst. Corrected!')
end
bin_width = 4; % in degrees
phase_avg_vid = zeros(size(adjvid,1), size(adjvid,2), 360/bin_width);
for ii = 1 : 360/bin_width
    inds = find(frame_lst(:, 4) >= -180+(ii-1)*bin_width & ...
                frame_lst(:, 4) < -180+ii*bin_width+4);
    inds(inds >= size(FoverFbar, 3)) = [];
    phase_avg_vid(:, :, ii) = mean(FoverFbar(:, :, inds), 3);
end

if GIFflag == 1
    tic
    making_GIF_file(PathName, [FileName(1:end-4), '-PhaseAvg'], phase_avg_vid)
    disp(['GIFation(!) took: ', num2str(toc)]);
end

%% Fourier Transform whole video
Fs = 1/mean(diff(frame_lst(:, 2)));  % Sampling frequency                    
L = length(frame_lst(:, 1));         % Length of signal
t = frame_lst(:, 2); % Time vector

Y = fft(FoverFbar, [], 3);
f = Fs*(0:(L-1))/L;

%%% finding stimulus frequency
stim_freq = find_stim_freq(t, frame_lst(:,4));
%stim_ind = find(round(f,3)==round(stim_freq,3));
[val, stim_ind] = min(abs(f - stim_freq));

Y_stim_real = 2*abs(Y(:, :, stim_ind)/L);
Y_stim_im = rad2deg(angle(Y(:, :, stim_ind)));

if rem(sweep_direction, 180) == 0 %%%transforming angles
    phase = cp_azdeg + Y_stim_im/pi*atan(SizeX/2/ScreenDist)*(-1)^(sweep_direction/180);
elseif rem(sweep_direction, 180) == 90
    phase = cp_eldeg + Y_stim_im/pi*atan(SizeY/2/ScreenDist)*(-1)^((sweep_direction-90)/180);
end

%% Plotting 
thr = 0.6; % Percentage of maximum amplitude
smoothed_amps = medfilt2(Y_stim_real);
smoothed_phase = wiener2(phase, [3,3]);
[r, c] = find(smoothed_amps>thr*max(smoothed_amps(:)));
colpixs = diag(smoothed_phase(r,c));

figure('Name', 'Real at stimulus freq', 'Position', [48,521,540,420])
imagesc(Y_stim_real); colorbar; hold on
imcontour(smoothed_amps, [thr*max(Y_stim_real(:)),thr*max(Y_stim_real(:))], 'k')
print([PathName, FileName(1:end-4), '-R@stimfreq.png'], '-dpng')

figure('Name', 'Imag at stimulus freq', 'Position', [48,521,540,420])
imagesc(phase); hold on; colormap('hsv')
caxis([min(colpixs), max(colpixs)]);colorbar
imcontour(smoothed_amps, [thr*max(Y_stim_real(:)),thr*max(Y_stim_real(:))], 'k')
print([PathName, FileName(1:end-4), '-I@stimfreq.png'], '-dpng')

%%%% phase with amplitude incorporated
figure('Name', 'Imag at stimulus freq', 'Position', [48,521,540,420])
img = imagesc(phase); hold on; colormap('hsv')
caxis([min(colpixs), max(colpixs)]);colorbar
namp = (Y_stim_real-min(Y_stim_real(:)))/(max(Y_stim_real(:))-min(Y_stim_real(:)));
img.AlphaData = namp;
imcontour(smoothed_amps, [thr*max(Y_stim_real(:)),thr*max(Y_stim_real(:))], 'k')
print([PathName, FileName(1:end-4), '-In@stimfreq.png'], '-dpng')

figure('Name', 'Imag at stimulus freq', 'Position', [48,521,540,420])
imagesc(phase); hold on; 
[Angles, col_map] = get_cmap(sweep_direction);
colormap(col_map)
if rem(sweep_direction, 180) == 90
    caxis([cp_eldeg-atan(SizeY/2/ScreenDist)*180/pi, cp_eldeg+atan(SizeY/2/ScreenDist)*180/pi]);colorbar
elseif rem(sweep_direction, 180) == 0
    caxis([cp_azdeg-atan(SizeX/2/ScreenDist)*180/pi, cp_azdeg+atan(SizeX/2/ScreenDist)*180/pi]);colorbar
end
imcontour(smoothed_amps, [thr*max(Y_stim_real(:)),thr*max(Y_stim_real(:))], 'k')
print([PathName, FileName(1:end-4), '-I2@stimfreq.png'], '-dpng')

%% Saving the workspace
supix = double(squeeze(mean(mean(adjvid(:,:,:)))));
save([PathName, FileName(1:end-4), '-Workspace.mat'], 'frame_lst', 'f', 'stim_freq','Y_stim_real', 'Y_stim_im', 'smoothed_amps', 'smoothed_phase', 'supix', 'phase')
