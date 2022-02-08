function c = calcontour(F1Contour,SIZ)

ncell = size(F1Contour,2);
for ii=1:ncell
    ccontour= contourc(reshape(F1Contour(:,ii),SIZ(1),SIZ(2)),[.05, .05]);
    ncontours = [find(ccontour(1,:)<1) size(ccontour,2)+1];
    [~,bigcontour] = max(diff(ncontours ));
    ccontour= ccontour(:,ncontours(bigcontour)+1:ncontours(bigcontour+1)-1);
    c.perim(ii) = sum(sqrt(diff(ccontour(1,2:end)).^2+diff(ccontour(2,2:end)).^2));
end
F1Contour(F1Contour<0.05)=0;
c.area = sum(F1Contour>0);
c.circ = c.perim.^2./(4*pi*c.area); % if elliptical c.circ~ (a/b+b/a)/2
