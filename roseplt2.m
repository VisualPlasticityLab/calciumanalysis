function h=roseplt2(figH,data1,data2);
global info;
THETA1=data1/info.steps(1)*2*pi;
THETA2=data2/info.steps(1)*2*pi;
steps=2*pi*[1/info.steps(1):1/info.steps(1):1];

h1=rose(figH,THETA1,steps);
hold on
h2=rose(figH,THETA2,steps);
h1.Color=[ 0.5 0 0];
h2.Color=[0 0 0.5]; 
legend('running','still','Orientation','horizontal');%'Location','northoutside',
