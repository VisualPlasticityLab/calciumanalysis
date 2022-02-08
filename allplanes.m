clear all
planes = dir('Plane?');
%%
p1f = dir([planes(1).name '\selected*\*_all.mat']);
p1f_ms = dir([planes(1).name '\selected*\*_as+ms*.mat']);
load(fullfile(p1f.folder,p1f.name),'F1Contour','F1cell','RGB1','eye1','eye2','fig1img');
load(fullfile(p1f_ms.folder,p1f_ms.name),'pair0','pair1','pair2');
cellnum = numel(pair0);

%% concetanate other planes
% 'PEAK0','PEAK00','PEAK2','ORI2','ODI0','ODI1',
for ii=2:numel(planes)
        p1f = dir([planes(ii).name '\selected*\*_all.mat']);
        p1f_ms = dir([planes(ii).name '\selected*\*_as+ms*.mat']);
        tmpf = load(fullfile(p1f.folder,p1f.name));
        tmpf_ms = load(fullfile(p1f_ms.folder,p1f_ms.name));
        F1Contour = cat(2,F1Contour,tmpf.F1Contour);
        F1cell = cat(2,F1cell,tmpf.F1cell);
        RGB1 = cat(2,RGB1,tmpf.RGB1);
        eye1 = catstr(eye1,tmpf.eye1);
        eye2 = catstr(eye2,tmpf.eye2);   
        fig1img(:,:,ii) = tmpf.fig1img(:,:,2);
        cellnum = [cellnum numel(tmpf_ms.pair0)];
        pair0 = cat(2,pair0,tmpf_ms.pair0);
        pair1 = cat(2,pair1,tmpf_ms.pair1);
        pair2 = cat(2,pair2,tmpf_ms.pair2);
end  
fig1img = fig1img/max(reshape(fig1img(200:end-200,200:end-200,:),[],1));
%%
h2=figure;imshow(fig1img);hold on;
SIZ = size(fig1img);
selected = find(pair0&pair1&pair2);
level = discretize(1:sum(cellnum),[1 cumsum(cellnum)+1]);
for jj= selected
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),[0.01 1],'LineColor',RGB1(:,jj),'LineWidth',3)
    text(F1cell(jj).Position(1),F1cell(jj).Position(2),F1cell(jj).String,...
        'Color',[level(jj)==1 level(jj)==2 level(jj)==3],'FontSize',16)
end
ODIcolorbar;
%%
saveas(h2, ['allplanes_selected' num2str(numel(selected)) '.png']);
saveas(h2, ['allplanes_selected' num2str(numel(selected)) '.fig']);
save('allplanes.mat','F1Contour','F1cell',...
    'RGB1','eye1','eye2','fig1img','pair0','pair1','pair2','cellnum','-v7.3')