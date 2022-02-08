function [tsamp,vsmooth] = realTrajectory(cvx,cvy,dT,dChannel,plotOn)


m1 = dChannel==0;
m2 = dChannel==1;


M1_X = cumsum(cvx(m1));
M2_X = cumsum(cvx(m2));
M1_Y = cumsum(cvy(m1));
M2_Y = cumsum(cvy(m2));


%%% resample data into even time intervals
M1_T = dT(m1);
M2_T = dT(m2);
dt = 0.1;
Tmax = min(max(M1_T),max(M2_T));
Tmin =max(min(M1_T),min(M2_T));

tsamp = Tmin:dt:Tmax;
M1_Xsamp = interp1(M1_T,M1_X,tsamp);
M2_Xsamp = interp1(M2_T,M2_X,tsamp);
M1_Ysamp = interp1(M1_T,M1_Y,tsamp);
M2_Ysamp = interp1(M2_T,M2_Y,tsamp);

M1_dXsamp = diff(M1_Xsamp);
M2_dXsamp = diff(M2_Xsamp);
M1_dYsamp = diff(M1_Ysamp);
M2_dYsamp = diff(M2_Ysamp);

%%% calculate real-world position by integrating head orientation
%%% and forward/orthogonal movement
C=9200;  %%% circumference in ticks 
scale_factor = C/(pi*20);  %%% ball is 20 cm diameter, so this is ticks/cm
plotOn = 1;

theta =zeros(size(M1_dYsamp));
x=zeros(size(M1_dYsamp)); y=zeros(size(M1_dYsamp));

length(M2_dXsamp);
find(isnan(M2_dXsamp));

dTheta=-2*pi*0.5*(M1_dXsamp+M2_dXsamp)/C;
M2_dYsamp =-1*M2_dYsamp;  %%% to keep axes following righthand rule
M1_dYsamp =-1*M1_dYsamp;  %%% because mouse is in front on 2p scope
theta(1)=0+dTheta(1);
x(1)=0+M1_dYsamp(1);
y(1) = 0+M2_dYsamp(1);
for t=2:length(M1_dYsamp);
    x(t) = x(t-1) + M1_dYsamp(t-1)*cos(theta(t-1)) - M2_dYsamp(t-1)*sin(theta(t-1));
    y(t) = y(t-1) + M1_dYsamp(t-1)*sin(theta(t-1)) + M2_dYsamp(t-1)*cos(theta(t-1));
    theta(t) = theta(t-1)+dTheta(t);
end


%%% plot trajectory
if plotOn
    figure
    plot(x/scale_factor,y/scale_factor)
    axis equal
end

if plotOn
    figure
    plot(theta);
    title('theta')
    ylabel('radians')
end

%%% plot velocity as a function of time
v = sqrt(diff(x).^2 + diff(y).^2)/(scale_factor*dt);
vsmooth = conv(v,ones(1,11))/10;
vsmooth=vsmooth(6:length(vsmooth)-4);
tsamp = tsamp(1:length(tsamp)-1)+dt/2;
v_thresh = 0.03*max(vsmooth);

if plotOn
    figure
    plot(tsamp,vsmooth,'g')
    xlabel('secs');
    ylabel('cm/sec');
end

total_distance = sum(v)*dt/100;
disp(sprintf('ran %f meters in %f mins',total_distance,max(tsamp)/60));
disp(sprintf('stationary = %f mins; moving = %f mins',sum(vsmooth<v_thresh)*dt/60,sum(vsmooth>v_thresh)*dt/60));
if sum((v>v_thresh)>0)  %%% to avoid divide by 0
    disp(sprintf('average moving speed = %f cm/sec \n ',mean(vsmooth(vsmooth>v_thresh))));
end