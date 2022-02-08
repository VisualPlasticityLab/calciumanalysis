function bad=pick(f)
%%%pick bad cells based on their center position%%%
%hcell=figure_enlarge(hcell,2);
try
    data=loadjson([strtok(f,'.') '.jmesh']);
catch
    [fn,path]=uigetfile('.jmesh','Load jmesh data');
    data= loadjson(fullfile(path,fn));
end
circle=data.jmesh{:};
total=numel(circle);

for n=1:total
    Co(:,n)=circle{n}.bbox;
end
%%%%%%Define bad signals: edge and background
ratio=[Co(2,:)-Co(1,:)]./[Co(4,:)-Co(3,:)];
elliptical=find(ratio>3|ratio<1/3);
%edge=([Co(3,:)+Co(4,:)])/2;
%out=find(edge<40|edge>720);
out=[];
defaultbad=cellstr(num2str([elliptical,out]));

%     selected=inputdlg('input bad cell#','Exclusion',[4 80],defaultbad);
%     bad=unique(str2num(selected{1}));
bad = [elliptical,out];

%    S_df=S_df(Cor);

