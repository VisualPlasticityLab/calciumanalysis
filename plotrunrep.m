function plotrunrep(matrix1,matrix2)
condition = {'day0','day1'};
numrep = size(matrix1,1);
nSteps = size(matrix1,2);
[names,names_dir] = num2ori(nSteps);

figure('Position',[600 800 800 300]);
subplot(1,2,1)
b1=bar(cat(1,numrep-sum(matrix1),sum(matrix1))','stacked');
b1(2).FaceColor='r';
ylabel('# of repetition')
xlabel('direction')
set(gca,'XTickLabel',names_dir);
%legend('still','run')
%legend('Boxoff')
title(condition{1},'FontSize',16)

subplot(1,2,2)
b2=bar(cat(1,numrep-sum(matrix2),sum(matrix2))','stacked');
%  legend({'still','run'},'orientation','Horizontal')
b2(2).FaceColor='r';
set(gca,'XTickLabel',names_dir);
% legend('Boxoff')
xlabel('direction')

title(condition{2},'FontSize',16)
saveas(gcf,'runplot.fig')
saveas(gcf,'runplot.png')

