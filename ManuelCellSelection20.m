clear all;close all;clc
%%
[day1sig,p1]=uigetfile('.signals','select signal file for day1eye');
s1 = matfile(fullfile(p1,day1sig));
fig1img =s1.V;
SIZ = size(fig1img);
F1Contour = s1.A;
sig = s1.sig;
%%
pos = strfind(strtok(day1sig,'c'),'_');
try
    load(fullfile(path1a,[figname1(1:pos(end)) 'memmap.mat']),'magnification');
    if mod(magnification,1)>0 || magnification ==2
        magnification = magnification*2.5;
    end
    load magnificationlist.mat ;
    mag = maglist(magnification);
catch
    mag = 2;
end
%%
h22=figure('Position',[1700 800 600 500]);
subplot('Position',[0.05 0.2 .9 .7])
imagesc(fig1img);axis off
colormap gray;
[pair0,pair1,c1]=autoselect(full(F1Contour),mag,sig(:,1:end-1));
selected1 = find(pair0&pair1);
caanalysisplot('.',selected1)
% Analysis Summary for automatic selection
% plotODIanalysis(fign,selected1)
% PLOTING individual TC, only draw potentially good cells within the selection region
%% Manuel inspection based on TC
pair2=zeros(size(pair1));
h1 = figure('Position',[0 800 600 400]);
h2 = figure('Position',[700 400 1600 200]);
for jj=selected1
    figure(h2); plot(sig(:,jj));
    axis off;title(sprintf('skewness=%.1f',skewness(sig(:,jj))));
    figure(h22);hold on;
    [~,hc]=contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),...
         [0.01 1],'LineColor',[1 0 0],'LineWidth',1)
     
    imgtitle = sprintf('reorganized_ page#%d.png',ceil(jj/1));
        figure(h1); imshow(imread((imgtitle)),'InitialMagnification',60);
        %set(gcf,'Position',[0 800 600 400]);
        choice = 1/menu('Good cell?','Y','Not gd','Noise');
        pair2(jj)= (choice==1);
        check(jj) = (choice ==.5);

    if pair2(jj) ==1
     hc.LineColor=[1 1 1];
    elseif check(jj) ==1
     hc.LineColor=[0 0 0];
    else
        hc.LineColor = 'none';
    end
end
%%
save('ManuelSelected','pair0','pair1','pair2','check','fig1img','F1Contour','sig','-v7.3');

%%
caanalysisplot('.',find(pair2))






%%
selectedtmp = find(pair0&pair1&pair2);
figure;imagesc(fig1img);colormap gray;hold on;
for jj=selectedtmp
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],...
        'LineColor','w','LineWidth',3)
end
for jj=selectedtmp
    hh2(jj)=text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
        'Color','k','FontSize',16);
end

plotODIanalysis(fign,selectedtmp)

%%
selected2 = find(pair0&pair1&pair2);
numcell = num2str(numel(selected2));
newfolder = ['selected_' numcell 'cells'];
mkdir(newfolder);
copyfile([fign 'all.mat'],newfolder);
cd(newfolder);

save([fign 'as+ms' numcell],'maskpic0','pair0','pair1','pair2','check','fign','numcell');
h2=figure;imshow(fig1img);hold on;
for jj=selected2
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
end
ODIcolorbar
saveas(h2,[fign 'as+ms' numcell '.png']);
plotODIanalysis(fign,selected2)


%% (Optional 1.) Manuel selection part of the imaging field 

figure(h2)
maskpic1=roipoly;
pair3 = reshape(maskpic0&maskpic1,1,prod(SIZ))*F1Contour >0; % the right half of the imaging field
save([fign 'as+ms' numcell],'pair3','maskpic1','-append');

selected3 = find(pair0&pair1&pair2&pair3);
h3=figure;imshow(fig1img);hold on;
for jj=selected3
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
end
ODIcolorbar
saveas(h3,[fign 'as+ms' num2str(numel(selected3)) '_RightPart.png']);
plotODIanalysis(fign,selected3);

selected4 = find(pair0&pair1&pair2&~pair3);
h4=figure;imshow(fig1img);hold on;
for jj=selected4
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
end
ODIcolorbar
saveas(h4,[fign 'as+ms' num2str(numel(selected4)) '_LeftPart.png']);
plotODIanalysis(fign,selected4);
%% (Optional 2.) Manuel inspection based on Red image
greenornot = F1Contour'*reshape(fig1img(:,:,2),[],1);
redornot = F1Contour'*reshape(fig1img(:,:,1),[],1)./greenornot;
figure;histogram(redornot,0:.1:1.5);
redmarkerlimit = prctile(redornot,75);
pair3 = pair0&pair1&pair2&(redornot'>redmarkerlimit);

%plot redcells in red and label in red
figure;imshow(fig1img(:,:,1));hold on;
for jj=selected2
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1-pair3(jj) 1-pair3(jj)],'LineWidth',3);
    text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
        'Color',[pair3(jj) 0 0],'FontSize',16)
end
%%
save([fign 'as+ms' numcell],'pair3','-append');
h3=figure;imshow(fig1img);hold on;
for jj=selected2
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1-pair3(jj) 1-pair3(jj)],'LineWidth',3);
end
saveas(h3,[fign 'as+ms' numcell '_redlabel.png']);
close(h3);

selected3 = find(pair0&pair1&pair2&pair3);
plotODIanalysis(fign,selected3);

selected4 = find(pair0&pair1&pair2&~pair3);
plotODIanalysis(fign,selected4);