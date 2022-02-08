function [pair0,pair1,pair2,c]=autoselect(img,F1Contour,mag,sig)
%% get pair0, pari1,pair2
[f,p]=uigetfile('*.mat','Load region from previous file?');
if f~= 0
    load(fullfile(p,f),'maskpic0');
    %load(fullfile(p,f),'pair0');
else
    maskpic0=roipoly;
    save('selectedregion.mat','maskpic0')
end
SIZ = size(img);
pair0 = reshape(maskpic0,1,[])*F1Contour >0;

% Set potential good cells criteria: response and size
c = calcontour(F1Contour,SIZ(1:2));
c.area_cutoff = mag^2*30 *[.55 3.5] ; % cellsize = [96 396] under 2X
c.circ_cutoff = 1.75; % 1/3 < a/b < 3 c.circ=5/3
c.skew_cutoff = 1;
if exist('sig','var')
    c.skew = skewness(sig);
    pair1 = c.area<=c.area_cutoff(2) & c.area>=c.area_cutoff(1)...
        & c.circ<=c.circ_cutoff & c.skew>c.skew_cutoff;
else
    pair1 = c.area<=c.area_cutoff(2) & c.area>=c.area_cutoff(1)...
        & c.circ<=c.circ_cutoff;
end

[~,values]= AutoCellSelection;%values=[corV;pV;peak]'
% pair22(:)=pair22(:)+1;
%% manuel inspection
h0 = figure();
imagesc(img,prctile(img(:),[1 99]));colormap gray;
hold on;
pair2(:) = zeros(size(pair22));
for jj=find(pair0&pair1&values(:,1)'>.1)
    figure(h0);hold on;
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),...
        [0.01 1],'LineColor',[1 0 0],'LineWidth',1)
    imgtitle = ['reorganized_ page#' num2str(jj) '.png'];
    h1=figure();imshow(imread([imgtitle]));
    set(h1,'Position',[0 200 1400 500])
    choice = menu(sprintf('cor=%.2f,P=%.2f,A=%.2f,Good?',values(jj,1),values(jj,2),values(jj,3)),'Y','N');
    pair2(jj)= (choice==1);
    close(h1);
    figure(h0);hold on;
    contour(1:SIZ(2),1:SIZ(1),reshape(F1Contour(:,jj),SIZ(1),SIZ(2)),...
        [0.01 1],'LineColor',pair2(jj)*[1 1 1],'LineWidth',3)
end
