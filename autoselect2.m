clear all
%%
[pair22,values]= AutoCellSelection;%values=[corV;pV;peak]'
pair2=zeros(size(pair22));
%%
% caanalysisplot('.',find(values(:,1)>.15&values(:,2)<.2),0)
selected= find(values(:,2)-values(:,1)<.05);
caanalysisplot('.',selected,0)
%% manuel inspection
for jj=selected(31:end)'
    imgtitle = ['reorganized_ page#' num2str(jj) '.png'];
    h1=figure();imshow(imread([imgtitle]));
    set(h1,'Position',[0 200 1400 500])
     prompt = sprintf('cor=%.2f,P=%.2f,A=%.2f,Good?',values(jj,1),values(jj,2),values(jj,3))
%     prompt = 'good?'
    choice = menu(prompt,'Y','N');
    pair2(jj)= (choice==1);
%     str = input(prompt,'s');
%     if strcmp(str,'y')
%         pair2(jj)=1;
%     end
    close(h1);
end

%%
figure;hold on
scatter(values(selected,1),values(selected,2),'bo');
scatter(values(pair2==1,1),values(pair2==1,2),'ro');
