function pieplt2(data1,data2);
hold on
h1 = pie(histcounts(data1));
(data1);
h2 = histcounts(data2);
h1.FaceColor=[ 0.5 0 0];
h2.FaceColor=[0 0 0.5]; 
legend('running','still','Location','north','Orientation','horizontal')