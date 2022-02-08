function fig=plotRG2(d1,matchedpair)

maxpicG = prctile(d1.ops.meanImg(:),[50 99]);
mimgG(:,:,1)=d1.ops.meanImg_chan2_corrected2/diff(maxpicG)*3;
mimgG(:,:,2)=(d1.ops.meanImg-maxpicG(1))/diff(maxpicG);
mimgG(:,:,3)=0;
cellct = 0;
bound =[];

fig=figure;
imshow(mimgG);
for M = matchedpair
    k = find(d1.iscell(:,1)==1,M);
    k = k(end);
    h = boundary(double(d1.stat{k}.xpix)',double(d1.stat{k}.ypix)');
    d1.stat{k}.xbound = d1.stat{k}.xpix(h);
    d1.stat{k}.ybound = d1.stat{k}.ypix(h);
    cellct = cellct + 1;
    bound{cellct} = [d1.stat{k}.xbound(:),d1.stat{k}.ybound(:)];
    figure(fig);hold on;

        if d1.redcell2(k,4)==1
        plot(bound{cellct}(:,1),bound{cellct}(:,2),'LineWidth',2,'Color',[1 0 0]);
        else
        plot(bound{cellct}(:,1),bound{cellct}(:,2),'LineWidth',1,'Color',[1 1 1]);
        end
        
    text(mean(bound{cellct}(:,1)),mean(bound{cellct}(:,2)),num2str(M),'color','k');
end