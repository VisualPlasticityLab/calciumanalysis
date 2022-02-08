clear all;close all;clc

tmpfile=dir('*sp*all.mat');
load(tmpfile.name);
% Analysis Summary for automatic selection
% plotODIanalysis(fign,selected1)
% PLOTING individual TC, only draw potentially good cells within the selection region
% plotTC(fign,selected1);
pair2=zeros(size(pair1));
check=zeros(size(pair1));
%% plotTC on ones with gd Fano factor 

%(done already)
% plotODIanalysis(fign,find(pair0&pair1));close all
plotTC(fign,find(pair0&pair1));close all
% caanalysisplot('..',1:numel(pair0),0);

% plotTC again only on neurons with gd max reponses
%plotODIanalysis(fign,find(pair0&pair1&(sum(PEAK1)>prctile(sum(PEAK1),50))));
%% Manuel inspection based on TC
RESTARTING = 1 %%%CAUSTION: UDPATE TO MAKE sure if you want to continue or starting from beginning
if RESTARTING ~= 1
    start_Middle = menu...
    ('Are you sure you want to start from the middle?',...
    'Y','No, start from the beginning');
    if start_Middle ~= 1
        RESTARTING = 1 
    end
end

h22=figure;
imshow(fig1img);
hold on;
for jj=find(pair0&pair1)
    
        if jj>=RESTARTING
    imgtitle1 = ['../reorganized_ page#' num2str(jj) '.png'];
    imgtitle = ['alltraces_cell#' num2str(jj) '_sp.png'];
    h0 = figure('Position',[100 200 150 850]);
    imshow(imread(imgtitle1));
    h1=figure('Position',[100 100 1300 400]);
    imshow(imread([imgtitle]));
    figure(h22);hold on;
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),...
        [0.01 1],'LineColor',[1 0 0],'LineWidth',1)
    try
        %         choice = menu(sprintf('cor=%.2f,P=%.2f,Good?',corV(jj),pV(jj)),'Y','N','check!');
        choice = menu('Good?','Y','N','check!');
        pair2(jj)= (choice==1);
        check(jj) = (choice ==3);
    catch
        pair2(jj) = 0;
    end
    figure(h22);hold on;
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),...
        [0.01 1],'LineColor',pair2(jj)*[1 1 1],'LineWidth',3)
    close(h1);close(h0);
        end
end

% %% automated grouping based on correlation
% [~,values]= AutoCellSelection;
% corV = values(:,1);
% pV = values(:,2);
% peak = values(:,3);
% try
%     figure;hold on;
%     scatter(corV(pair2==0),pV(pair2==0),'ro');
%     scatter(corV(pair2==1),pV(pair2==1),'bo');
%     legend('bad cells','good cells')
%     figure;
%     histplt2(corV(pair2==0),corV(pair2==1),-.025:.05:1.025);
%     legend({'bad cells','good cells'})
% end
% %% check
% selectedtmp = find(pair0&pair1&pair2);
% % plotODIanalysis(fign,selectedtmp)
% 
% figure;imshow(fig1img);colormap gray;hold on;
% for jj=selectedtmp
%     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
% end
% for jj=selectedtmp
%     hh2(jj)=text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
%         'Color',RGB1(:,jj),'Color','k','FontSize',16);
% end
% ODIcolorbar

%% SAVE
selected2 = find(pair0&pair1&pair2);
numcell = num2str(numel(selected2));
newfolder = ['selected_' numcell 'cells'];
mkdir(newfolder);
copyfile([fign 'all.mat'],newfolder);
cd(newfolder);
if ~exist('check','var')
    check =[];
end
save([fign 'as+ms' numcell],'maskpic0','pair0','pair1','pair2','check','fign','numcell');
plotODIanalysis(fign,selected2)

h2=figure;imshow(fig1img);hold on;
for jj=selected2
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
end
ODIcolorbar
saveas(h2,[fign 'as+ms' numcell '.png']);



%% skip for now,JS -2020 May

% %% (Optional 1.) Manuel selection part of the imaging field
% 
% figure(h2)
% maskpic1=roipoly;
% pair3 = reshape(maskpic0&maskpic1,1,prod(SIZ))*F1Contour >0; % the right half of the imaging field
% save([fign 'as+ms' numcell],'pair3','maskpic1','-append');
% 
% selected3 = find(pair0&pair1&pair2&pair3);
% h3=figure;imshow(fig1img);hold on;
% for jj=selected3
%     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
% end
% ODIcolorbar
% saveas(h3,[fign 'as+ms' num2str(numel(selected3)) '_RightPart.png']);
% plotODIanalysis(fign,selected3);
% 
% selected4 = find(pair0&pair1&pair2&~pair3);
% h4=figure;imshow(fig1img);hold on;
% for jj=selected4
%     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
% end
% ODIcolorbar
% saveas(h4,[fign 'as+ms' num2str(numel(selected4)) '_LeftPart.png']);
% plotODIanalysis(fign,selected4);
% %% (Optional 2.) Manuel inspection based on Red image
% greenornot = F1Contour'*reshape(fig1img(:,:,2),[],1);
% redornot = F1Contour'*reshape(fig1img(:,:,1),[],1)./greenornot;
% figure;histogram(redornot,0:.1:1.5);
% redmarkerlimit = prctile(redornot,75);
% pair3 = pair0&pair1&pair2&(redornot'>redmarkerlimit);
% 
% %plot redcells in red and label in red
% figure;imshow(fig1img(:,:,1));hold on;
% for jj=selected2
%     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1-pair3(jj) 1-pair3(jj)],'LineWidth',3);
%     text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
%         'Color',[pair3(jj) 0 0],'FontSize',16)
% end
% %%
% save([fign 'as+ms' numcell],'pair3','-append');
% h3=figure;imshow(fig1img);hold on;
% for jj=selected2
%     contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',[1 1-pair3(jj) 1-pair3(jj)],'LineWidth',3);
% end
% saveas(h3,[fign 'as+ms' numcell '_redlabel.png']);
% close(h3);
% 
% selected3 = find(pair0&pair1&pair2&pair3);
% plotODIanalysis(fign,selected3);
% 
% selected4 = find(pair0&pair1&pair2&~pair3);
% plotODIanalysis(fign,selected4);