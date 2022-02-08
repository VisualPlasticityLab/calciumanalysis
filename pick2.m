function [bad,h]=pick2(sigT,win_sig,Cor1)

[~,peak,~]=cal_ER(sigT,win_sig);%peak:1*Var*ncell
peak=squeeze(peak);
Var=size(peak,1);
ncell=size(peak,2);


h=figure('Name','all peak values for each cell','position',[000 200 2000 600]);
h0=subplot(2,1,1);hold on;
x=ones(Var,1)*Cor1;
scatter(x(:),peak(:),'o','jitter','on', 'jitterAmount', 0.2);
% 
nnz=sum(peak>0);
ym=mean(peak);ym(isnan(ym))=0;

peak(peak==0)=NaN;
ynz=nanmean(peak);ynz(isnan(ynz))=0;
%ymax=prctile(peak,90);ymax(isnan(ymax))=0;
% 
% 
h0.XTick=Cor1;
%yr=std(squeeze(peak));
plot(Cor1,ym(:),'*');
plot(Cor1,ynz(:),'x');
xlabel(' cell #');
ylabel('max dF/F for each stim');
% 
% h1=subplot(2,3,4);hold on;H1=histogram(nnz,20);xlabel('# of nonzeros')
% h2=subplot(2,3,5);hold on;H2=histogram(ynz,20);xlabel('mean of nonzero values')
% h3=subplot(2,3,6);hold on;H3=histogram(ym,40);xlabel('mean of all values')
% 
% 
% % prompt = {'Enter cutoff bins'};
% % dlg_title = 'Parsing';
% % num_lines = [1 80];
% % defaultans = {'1 1 1'};
% % answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
% % Edge=str2num(answer{1});
% Edge=[1 1 1 ];
% nnz_thr=H1.BinEdges(Edge(1)+1);
% mnz_thr=H2.BinEdges(Edge(2)+1);
% mean_thr=H3.BinEdges(Edge(3)+1);
% plot(h0,[ 0 Cor1(end)],ones(2,1)*[mnz_thr mean_thr],'--');
% plot(h1,[nnz_thr nnz_thr],[0 max(H1.Values)],'r--')
% plot(h2,[mnz_thr mnz_thr],[0 max(H2.Values)],'r--')
% plot(h3,[mean_thr mean_thr],[0 max(H3.Values)],'r--')


%out=find(nnz<=nnz_thr&ynz<mnz_thr&ym<mean_thr);
% out=find(nnz<=nnz_thr);
% defaultbad=cellstr(num2str(Cor1(out)));
% selected=inputdlg('input bad cell#','Exclusion',[4 80],defaultbad);
% bad=unique(str2num(selected{1}));
bad = find(nanmax(peak,[],1)<.1);
disp(['Found bad cells: ' num2str(bad) ])
% peak:1*Var*ncell
    