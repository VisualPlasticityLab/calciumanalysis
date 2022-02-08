%% LOAD green and red images
% if you have only one plane change the steps from 4 to 2
close all;clear all
plane_num=0 ; %plane number. either 0 (plane1) , 1 (plane2),2(plane 3)
nplane =3;
[fname,fpath] = uigetfile('*.tif');

ntiff = numel(imfinfo([fpath, fname]));
fprintf('Loading Green&Red...\n')
samplingrate = 10; 
green1_img = zeros(512, 796, floor(ntiff/nplane/2/samplingrate));
red1_img = zeros(512, 796, floor(ntiff/nplane/2/samplingrate));
for ii = 1:size(green1_img,3) % 1=green plane 1, 2=red plane 1, 3=green plane 2 ,4=red plane2
    green1_img(:, :, ii) = imread([fpath, fname], (ii-1)*2*nplane*samplingrate+plane_num+1);
    red1_img(:, :, ii) = imread([fpath, fname],  (ii-1)*2*nplane*samplingrate+plane_num+2);
end

%%% LOAD suite2p output
load(['suite2p/plane', num2str(plane_num), '/Fall.mat'],'stat');

%% START analyzing
r_1 = 15;
r_2 = 10;
avg_img = mean(green1_img, 3);

[X, Y] = meshgrid(1:796, 1:512);

gIntensity = [];
rIntensity = [];
centers = [];
nFrames = size(green1_img, 3);
cell_id = [];
cnter = 0;
for cid = 1 : length(stat)
%     if iscell(cid, 1) == 1
        cnter = cnter + 1;
        h = boundary(double(stat{cid}.xpix)',double(stat{cid}.ypix)');
        xbound = stat{cid}.xpix(h);
        ybound = stat{cid}.ypix(h);
        xcenter = mean(xbound);
        ycenter = mean(ybound);
        
        distm = ((X-xcenter).^2 + (Y-ycenter).^2).^0.5;
        distm(distm <= r_1) = 1;
        distm(distm > r_1) = nan;
        
         aux = avg_img.*distm;
        [ycenter, xcenter] = find(aux == max(aux(:)));
        distm = ((X-xcenter).^2 + (Y-ycenter).^2).^0.5;
        distm(distm <= r_2) = 1;
        distm(distm > r_2) = nan;
        
%         figure(101); imagesc(avg_img.*distm);
        %choice = menu('Include this cell?','YES','NO');
%         if iscell(cid, 1) == 1
            centers = [centers; xcenter, ycenter];
            gvec = nan(nFrames, 1);
            rvec = nan(nFrames, 1);
            for frame = 1 : nFrames
                gvec(frame) = nanmean(nanmean(squeeze(green1_img(:, :, frame)).*distm));
                rvec(frame) = nanmean(nanmean(squeeze(red1_img(:, :, frame)).*distm));
            end
             cell_id = [cell_id, cnter];
             gIntensity = [gIntensity, gvec];
             rIntensity = [rIntensity, rvec];
        end
%     end
% end
 
%close 
%% plotting all cell traces
% for ii = 1 : size(rIntensity, 2)
%     figure(202); clf
%     plot(gIntensity(:, ii), 'g'); hold on
%     plot(rIntensity(:, ii), 'r');
%     title(['cell id: ', num2str(cell_id(ii))])
%     pause
% end


% 
% %% plotting individual cell traces
% cid = 11;
% for ii = find(cell_id==cid)
%     figure(111); clf
%     plot(gIntensity(:, ii), 'g'); hold on
%     plot(rIntensity(:, ii), 'r');
%     title(['cell id: ', num2str(cell_id(ii))])
% end

%% plotting all cells scatter plots with fit
% for ii = 1 : size(rIntensity, 2)
%     figure (404);clf
%   
%     plot(gIntensity(:, ii), rIntensity(:, ii), 'o', 'MarkerEdgeColor', 'k', ...
%                 'MarkerFaceColor', 'k', 'MarkerSize', 3); alpha 0.4; hold on
%     
%     [fitres, gof] = createFit(gIntensity(:, ii), rIntensity(:, ii));
%     
%     plot([min(gIntensity(:,ii)), max(gIntensity(:,ii))], ...
%          fitres.p2 + fitres.p1*([min(gIntensity(:,ii)),max(gIntensity(:,ii))]),'r-', 'linewidth', 2);
%     title(['cell id:', num2str(cell_id(ii)), ' , R-squ: ', num2str(gof.rsquare)]); box off
%     xlabel('Green channel')
%     ylabel('Red channel')
%     pause
% end
%%
% 
   %% plotting individual scatter plots with fits
% cid = 22;
% for ii = find(cell_id==cid)
%     figure (112);clf
%   
%     plot(gIntensity(:, ii), rIntensity(:, ii), 'o', 'MarkerEdgeColor', 'k', ...
%                 'MarkerFaceColor', 'k', 'MarkerSize', 3); alpha 0.4; hold on
%     
%     [fitres, gof] = createFit(gIntensity(:, ii), rIntensity(:, ii));
%     
%     plot([min(gIntensity(:,ii)), max(gIntensity(:,ii))], ...
%          fitres.p2 + fitres.p1*([min(gIntensity(:,ii)),max(gIntensity(:,ii))]),'r-', 'linewidth', 2);
%     title(['cell id:', num2str(cell_id(ii)), ' , R-squ: ', num2str(gof.rsquare)]); box off
%     xlabel('Green channel')
%     ylabel('Red channel')
% end
%% 
 %% plotting individual scatter plots
% cid = 29;
% for ii = find(cell_id==cid)%1 : size(rIntensity, 2)
%    
%     figure(299); clf
%     %subplot(6, 6, ii)
%     plot(gIntensity(:, ii), rIntensity(:, ii), 'o');
%     xlabel('green'); ylabel('red')
%     title(['cell id: ', num2str(cell_id(ii))])
% end


%% plotting correctd red image
corrected_img = zeros(size(avg_img));
averages = nan(1, size(rIntensity, 2));
ratio = nan(1, size(rIntensity,2));
isred = nan(1, size(rIntensity,2));
fprintf('adjusting ceell intensity: ') 
for ii = 1 : size(rIntensity, 2)
    fprintf([num2str(ii), '\n'])
    fitres=fit(gIntensity(:, ii),rIntensity(:, ii),'poly1');
%    plot(f, [0; gIntensity(:, ii)],[f(0); rIntensity(:, ii)]);
      isred(:,ii) = fitres(0);
    distm = ((X-centers(ii, 1)).^2 + (Y-centers(ii,2)).^2).^0.5;
    distm(distm <= r_1) = 1;
    distm(distm > r_1) = 0;
    
    aux = mean(red1_img - fitres.p1*green1_img, 3); % adjusted red for each ROI
    R = sum(sum(aux.*distm))/mean(green1_img(:)); % adjusted red normalized for green brightness
    averages(ii) = sum(sum(aux.*distm))/pi/r_2^2; % average adjusted red intensity at each ROI
    ratio(ii) = R/pi/r_2^2;
    corrected_img = corrected_img + aux.*distm;
    
end
%%
figure();
colmap = [linspace(1,0,1e3)', ones(1000,1), linspace(1,0,1e3)'];
corrected_img(corrected_img==0) = nan;
ax1 = subplot(2,2,1);
p1 = imagesc(mean(green1_img,3)); colormap(colmap); colorbar
title('Avg green')
ax2 = subplot(2,2,2);
p2 = imagesc(mean(red1_img,3)); colormap(colmap); colorbar
title('Avg red')
ax3 = subplot(2,2,3);
p3 = imagesc(corrected_img); colormap(colmap); colorbar
title('Avg adjusted red')
for ii = 1 : size(centers, 1)
    text(centers(ii,1), centers(ii,2), num2str(cell_id(ii)), 'FontSize',6)%, ',', ...
%         num2str(round(averages(ii),0))], 'Color', 'k' , 'FontSize',5)
end 
ax4 = subplot(2,2,4);
p4 = imagesc(corrected_img/mean(green1_img(:))); colormap(colmap); colorbar
linkaxes([ax3,ax4], 'xy')
for ii = 1 : size(centers, 1)
    text(centers(ii,1), centers(ii,2), num2str(cell_id(ii)),'FontSize',6) %, ',', ...
%         num2str(round(ratio(ii),2))], 'Color', 'k', 'FontSize',5)
end
title('Avg adjusted red divid. by Avg green')

saveas(gcf,['green_red_corrected_plane',num2str(plane_num+1),'.fig'])


%% increase contrast
figure()
imagesc(corrected_img); colormap(jet)
for ii = 1 : size(centers, 1)
    text(centers(ii,1), centers(ii,2), num2str(cell_id(ii)), 'FontSize',6)%, ',', ...
%         num2str(round(averages(ii),0))], 'Color', 'k' , 'FontSize',5)
end 

saveas(gcf,['corrected_img_plane',num2str(plane_num+1),'.fig'])
saveas(gcf,['corrected_img_plane',num2str(plane_num+1),'.tif'])


%% make contrast image

I = imread(['corrected_img_plane',num2str(plane_num+1),'.tif']);
B = I(:,:,2);
figure()
J = adapthisteq(B,'clipLimit',0.02,'Distribution','rayleigh');
imshow(J)

saveas(gcf,['corrected_img_contrast_plane',num2str(plane_num+1),'.tif'])
% sort by ratio
xx=[cell_id',ratio'];[~,idx]=sort(xx(:,2));xx(idx,:)


%% increase contrast for responsive cells
% figure()
% imagesc(corrected_img); colormap(parula)
% resp_cells= [1
% 2
% 5
% 6
% 7
% 8
% 9
% 10
% 11
% 13
% 14
% 15
% 16
% 17
% 18
% 19
% 20
% 21
% 23
% 24
% 25
% 26
% 29
% 31
% 32
% 33]
% 
% for ii = 1 : length(resp_cells)
%     fprintf('ii: %d, resp: %d, c: %d, cc: %d \n', ii, resp_cells(ii), centers(ii,1), centers(resp_cells(ii),1))
%     text(centers(resp_cells(ii),1), centers(resp_cells(ii),2), num2str(cell_id(resp_cells(ii))), 'FontSize',14)%, ',', ...
% %         num2str(round(averages(ii),0))], 'Color', 'k' , 'FontSize',5)
% end 
% 
% saveas(gcf,['corrected_img_respcell_plane_',num2str(plane_num+1),'.fig'])
% saveas(gcf,['corrected_img_respcell_plane_',num2str(plane_num+1),'.tif'])

% 
% %% make contrast image
% 
% I = imread(['corrected_img_respcell_plane_',num2str(plane_num+1),'.tif']);
% B = I(:,:,2);
% figure()
% J = adapthisteq(B,'clipLimit',0.02,'Distribution','rayleigh');
% imshow(J)
% 
% 
% 
% 
% saveas(gcf,['corrected_img_contrast_respcell_plane',num2str(plane_num+1),'.tif'])
% % sort by ratio
% xx=[cell_id',ratio'];[~,idx]=sort(xx(:,2));xx(idx,:)
