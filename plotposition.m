%function plotposition

%%
if ~exist('day1peakf')
[day1peakf,path1]=uigetfile('peakSI.mat','select peakSI matfile for day1eye');
[day2peakf,path2]=uigetfile('peakSI.mat','select peakSI matfile for day2eye');
day1peak=load(fullfile(path1,day1peakf));
day2peak=load(fullfile(path2,day2peakf));

[alignedfig,figpath]=uigetfile('*.fig','find the aligned fig');
h4 = openfig(fullfile(figpath,alignedfig));
SIZ = size(getimage(h4));
end
if ~exist('F1Contour')
[contourfile,contourfilepath]=uigetfile('*&*_selected*.mat',...
    'find the mat file with cell contours');
load(fullfile(contourfilepath,contourfile));
end

[prefilename,path]=uigetfile('selected_*pairs.mat'...
    ,'Pick the pre-matched pairs?');
pre = load(fullfile(path,prefilename));
goodcell1 = pre.selectedpair(:,1);
goodcell2 = pre.selectedpair(:,2);
%%
plotrunrep(day1peak.matrix,day2peak.matrix)
%%
ori = (size(pre.peak_selected,2)-1)/2;
colortype = ori;
h = (1:colortype)'/colortype;
s = ones(colortype,1);
v = ones(colortype,1);
clut =squeeze(hsv2rgb(h,s,v));
clut(ori+1:ori*2,:) = clut/2;
clut(ori*2+1,:) = [1 1 1]/2;
%clut2 =clut/2;

de=0:180/ori:360;
dx = cosd(de)*10;
dy = sind(de)*10;
dx(end) = 0; dy(end) = 0; % noise level sets to 0

%% define the condition run/still/both to plot
cond = 'day2';
%%
switch cond
    case 'both'
        fit_selected = combinestatus_fitselected(pre.fit_selected,day1peak.matrix,day2peak.matrix);
        peak_selectedday = combinestatus(pre.peak_selected,day1peak.matrix,day2peak.matrix);
    case 'run'
        peak_selectedday = pre.peak_selected([2 4],:,:);
    case 'still'
        peak_selectedday = pre.peak_selected([1 3],:,:);  
    case 'day1'
        peak_selectedday = pre.peak_selected([1 2],:,:);
    case 'day2'
        peak_selectedday = pre.peak_selected([3 4],:,:);  
end

peak= discretize((max(peak_selectedday(1,:,:))),3);

x=[];y=[];
sz1 = peak(:);
sz2=[];
dir1=[];
dir2=[];
nCell = numel(goodcell1);

% get each cell's position
for j=1:nCell
    cell_1 = goodcell1(j);
    cell_2 = goodcell2(j);
    [x1,y1]=find(reshape(F1Contour(:,cell_1),SIZ(1),SIZ(2))==max(F1Contour(:,cell_1)),1);
    x=cat(1,x,x1);
    y=cat(1,y,y1);
end

% get size value
switch cond
    case 'both'
        for j=1:nCell
            sz2=cat(1,sz2,fit_selected{j}.a2/fit_selected{j}.a1*sz1(j));
            dir1(end+1) = fit_selected{j}.ori1;
            dir2(end+1) = fit_selected{j}.ori2;
            cond1='day1';
            cond2='day2';
        end
    case 'run'
        for j=1:nCell
            sz2=cat(1,sz2,fit_selected{j}.a2R/fit_selected{j}.a1R*sz1(j));
            dir1(end+1) = fit_selected{j}.ori1R;
            dir2(end+1) = fit_selected{j}.ori2R;
            cond1='day1';
            cond2='day2';
        end
    case 'still'
        for j=1:nCell
            sz2=cat(1,sz2,fit_selected{j}.a2S/fit_selected{j}.a1S*sz1(j));
            dir1(end+1) = fit_selected{j}.ori1S;
            dir2(end+1) = fit_selected{j}.ori2S;
            cond1='day1';
            cond2='day2';
        end
    case 'day1'
        for j=1:nCell
            sz2=cat(1,sz2,fit_selected{j}.a1R/fit_selected{j}.a1S*sz1(j));
            dir1(end+1) = fit_selected{j}.ori1S;
            dir2(end+1) = fit_selected{j}.ori1R;
            cond1='still';
            cond2='run';
        end
     case 'day2'
        for j=1:nCell
            sz2=cat(1,sz2,fit_selected{j}.a2R/fit_selected{j}.a2S*sz1(j));
            dir1(end+1) = fit_selected{j}.ori2S;
            dir2(end+1) = fit_selected{j}.ori2R;
            cond1='still';
            cond2='run';
        end
end

%     ct1 = text(F1cell(cell_1).Position(1),F1cell(cell_1).Position(2),F1cell(cell_1).String,'Color',clut(fit_selected{j}.ori1,:),'FontSize',8);
%     cc1 = contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,cell_1),SIZ(1),SIZ(2)),[0.01 1],'LineColor',clut(fit_selected{j}.ori1,:));
%     ct2 = text(F2cell(cell_2).Position(1),F2cell(cell_2).Position(2),F2cell(cell_2).String,'Color',clut(fit_selected{j}.ori2,:),'FontSize',8);
%     cc2 = contour(1:SIZ(2),1:SIZ(1),reshape(F2Contour(:,cell_2),SIZ(1),SIZ(2)),[0.01 1],'LineColor',clut2(fit_selected{j}.ori2,:));
% 
%      plot(x1,y1,'o','Color',clut(fit_selected{j}.ori2,:),'LineWidth',60*fit_selected{j}.a2,...
%          'Markersize',60*fit_selected{j}.a1,'MarkerFaceColor',clut(fit_selected{j}.ori1,:))
% 
%     plot(x1,y1,'o','Color',clut(fit_selected{j}.ori2,:),...
%         'Markersize',12*fit_selected{j}.a2/fit_selected{j}.a1,'MarkerFaceColor',clut(fit_selected{j}.ori2,:))

%1stday
%     plot(x1,y1,'o','Color',clut(fit_selected{j}.ori1,:),...
%         'Markersize',6*peak(j),'MarkerFaceColor',clut(fit_selected{j}.ori1,:))

h5 = figure('Position', [100 200  1000 600]);hold on
hh6 = gca;
% hh6.XLim =[0 850];
% hh6.YLim =[0 850];
hh1=scatter(x,y,200*sz1,clut(dir1,:),'filled') % day1 orientation
hh2=scatter(x,y,240*sz2,clut(dir2,:),'LineWidth',2)  % day2 orientation
alpha(hh1,.3);
% axis equal

axis off
%scalebar('ScaleLength',[36,64],'Unit',{'20mm','20mm'},'Location','southwest')
title(cond)
%
% h6 = figure('Position', [100 200  1500 1000]);hold on
for j=1:nCell
    arrow([x(j) y(j)],[x(j)+dx(dir1(j))*sqrt(sz1(j)) y(j)+dy(dir1(j))*sqrt(sz1(j))],20,...
        'EdgeColor',clut(dir1(j),:),'FaceColor',[1 1 1],'LineWidth',1,'Length',8,'tipangle',30);
    arrow([x(j) y(j)],[x(j)+dx(dir2(j))*sqrt(sz2(j)) y(j)+dy(dir2(j))*sqrt(sz2(j))],20,...
        'EdgeColor',clut(dir2(j),:),'FaceColor',clut(dir2(j),:),'LineWidth',1,'Length',8,'tipangle',30);
end


% axis tight
x0 = hh6.XLim(2);
y0 = hh6.YLim(2);
oricolormap(x0-10*ori,y0-10*ori,10*ori,ori,cond1,cond2)
% oricolormap(150,650,8*ori,ori);
alpha(.3)

% set(gca,'YDir','reverse')
saveas(h5,[cond '_postion.fig'])
saveas(h5,[cond '_postion.png'])
