d1 = load('Fall.mat');
d1.F_chan2= d1.Fcell_chan2;

nROI=size(d1.iscell,1);
if size(d1.redcell,1)<nROI
    d1.redcell(end+1:nROI,:)=0;
    d1.F_chan2(end+1:nROI,:)=0;
end

%%  linear regression is generated, mean Image corrected,NOTE stat.xpix and stat.ypix are reversed 
debug=0;
debug2=0;
d1.ops.meanImg_chan2_corrected2 = fitmin(d1.ops.meanImg,d1.ops.meanImg_chan2,0);
d1.ops.meanImg_chan2_corrected3 = d1.ops.meanImg_chan2_corrected2;

SIZ = size(d1.ops.meanImg_chan2);
Fsig_prcF = zeros(nROI,1);
Rsig_meanF = zeros(nROI,2);
Rsig_meanImg = zeros(nROI,3);
for i=1:nROI
%      figure(i)
    [d1.F_chan2_corrected(i,:),p] = fitmin(d1.F(i,:),d1.F_chan2(i,:),debug);
%     d1.Fneu_chan2_corrected(i,:) = fitmin(d1.Fneu(i,:),d1.Fneu_chan2(i,:),debug2);
%      title(sprintf('cell%d:red %d',i,round(mean(d1.F_chan2_corrected(i,:)))));
    
%     Fit_mean(i,:d1) = polyfit(d1.green(i,:),d1.red(i,:),1);
%     d1.red2(i,:)=d1.red(i,:)-polyval(Fit_min(i,:),d1.green(i,:));
    Fsig_prcF(i) = prctile(d1.F(i,:),95) ;
    Rsig_meanF(i,:) = [mean(d1.F_chan2_corrected(i,:)) mean(d1.F_chan2(i,:))];
    
    ROIidx = sub2ind(SIZ,d1.stat{i}.ypix,d1.stat{i}.xpix); %
    d1.ops.meanImg_chan2_corrected3(ROIidx) = d1.ops.meanImg_chan2(ROIidx) -polyval(p,d1.ops.meanImg_chan2(ROIidx));
    Rsig_meanImg(i,:) = [mean(d1.ops.meanImg_chan2_corrected3(ROIidx)) mean(d1.ops.meanImg_chan2_corrected2(ROIidx)) mean(d1.ops.meanImg_chan2_corrected(ROIidx))];
end

%% identify red cells
% print default analysis
cellidx = find(d1.iscell(:,1)==1)';
redcellidx=find(d1.redcell(:,1)==1)';
excellidx = setdiff(cellidx,redcellidx);
noncellidx = setdiff(1:nROI,cellidx);

figure;hold on;
scatter(Fsig_prcF(noncellidx),Rsig_meanImg(noncellidx,1),'bx');
scatter(Fsig_prcF(cellidx),Rsig_meanImg(cellidx,1),'ko');
scatter(Fsig_prcF(redcellidx),Rsig_meanImg(redcellidx,1),'ro');
xlabel('Green signal')
ylabel('Red signal from corrected red mean image')
legend({'non-cell','ex cell','red cell'})

% recalculate using our method
x = double(Fsig_prcF);
y = double(Rsig_meanImg(:,1));

[p,S] = polyfit(x(y<prctile(y,90)),y(y<prctile(y,90)),1); 
[y_fit,delta] = polyval(p,x,S);

%Plot the original data, linear fit, and 95% prediction interval y±2Δ.

plot(x,y_fit,'r-')
plot(x,y_fit+3*delta,'m--',x,y_fit-3*delta,'m--')
title('Linear Fit of Data with 99% Prediction Interval')
legend({'non-cell','ex cell','red cell','Linear Fit','99% Prediction Interval'})

isred = find(y>y_fit+3*delta)';%find cells that is likely to be red
newredcellidx = setdiff(intersect(isred,cellidx),redcellidx);

%% Plot image using default redchannel, default corrected redchannel, or our corrected redchannel

figure(130)
plotImgwCell(d1.ops.meanImg,d1.ops.meanImg_chan2_corrected,d1.stat,excellidx,redcellidx)
title('default corrected redchannel')

figure(133)
plotImgwCell(d1.ops.meanImg,d1.ops.meanImg_chan2_corrected2,d1.stat,excellidx,redcellidx,newredcellidx)
title('our corrected redchannel')

%%
figure(140)
plotImgwCell(ops_pre{1}.meanImg,ops_post{1}.meanImg_chan2_corrected,d1.stat,excellidx,redcellidx)
%%
figure(141)
plotImgwCell(ops_pre{1}.meanImg,ops_pre{1}.meanImg_chan2_corrected,d1.stat,excellidx,redcellidx)

%%
figure(129)
plotImgwCell(d1.ops.meanImg,d1.ops.meanImg_chan2_corrected)


figure(130)
plotImgwCell(ops_post{1}.meanImg,ops_post{1}.meanImg_chan2_corrected,d1.stat,excellidx,redcellidx)

figure(131)
plotImgwCell(ops_pre{1}.meanImg,ops_pre{1}.meanImg_chan2_corrected)

figure(132)
plotImgwCell(ops_post{1}.meanImg,ops_post{1}.meanImg_chan2_corrected2,stat_post,1:101,find(redcell(1:101,1)==1))


%%
boxplot([d1.ops.meanImg(:),ops_post{1}.meanImg(:),ops_pre{1}.meanImg(:)])
boxplot([d1.ops.meanImg_chan2(:),ops_post{1}.meanImg_chan2,ops_pre{1}.meanImg_chan2])
boxplot([d1.ops.meanImg_chan2_corrected(:),ops_post{1}.meanImg_chan2_corrected(:),ops_pre{1}.meanImg_chan2_corrected(:)])

     %%  
            prompt = {'delete old ones','add new cells'};
            figtitle = 'update red cells';
            cells=inputdlg(prompt, figtitle,[2 80]);
            rdcellidx = union(newredcellidx ,str2num(cells{1}));
            rdcellidx = setdiff(rdcellidx,str2num(cells{2}));


%% plot figure with red cells
function plotImgwCell(gMeanImg,rMeanImg,stat,cellidx,rdcellidx,newcellidx)
maxpicR = prctile(rMeanImg(:),99.9);
maxpicG = prctile(gMeanImg(:),99);
mimg(:,:,1)=rMeanImg/maxpicR;
mimg(:,:,2)=gMeanImg/maxpicG;
mimg(:,:,3)=0;
imshow(mimg,'InitialMagnification',200);
hold on;
if exist('cellidx','var')
for k = cellidx
    h = boundary(double(stat{k}.xpix)',double(stat{k}.ypix)');
    xbound = stat{k}.xpix(h);
    ybound = stat{k}.ypix(h);
    plot(xbound(:),ybound(:),'color','r','LineWidth',.5);   
    text(mean(xbound),mean(ybound),num2str(k),'color','k');
end
end

if exist('rdcellidx','var')
 
for k = rdcellidx
    h = boundary(double(stat{k}.xpix)',double(stat{k}.ypix)');
    xbound = stat{k}.xpix(h);
    ybound = stat{k}.ypix(h);
    plot(xbound(:),ybound(:),'color','r','LineWidth',2);
    text(mean(xbound),mean(ybound),num2str(k),'color','r');
end
end

if exist('newcellidx','var')
    for k = newcellidx
        h = boundary(double(stat{k}.xpix)',double(stat{k}.ypix)');
        xbound = stat{k}.xpix(h);
        ybound = stat{k}.ypix(h);
        plot(xbound(:),ybound(:),'color','r','LineWidth',.5);
        text(mean(xbound),mean(ybound),num2str(k),'color','w');
    end
end

end