clear;clc;
addpath('./CSV')
%% Define parameter

h = 0.005;

lf = 0.8415;
lr = 0.6885;
m = 270;
Iz = 105;
cf = 18748;
cr = 15339;

vx_threshold = 2.8;



%% Define Matrix
P = diag([10 10 10]);
% R = diag([0.5 0.5]);
% Q = diag([0 0 0]);

R = diag([1000000 10000]);
Q = diag([10 10 40]);

%This is good 
% R = diag([2000 500]);
% Q = diag([10 20 30]);

%This is for non-flip data
% R = diag([2000 500]);
% Q = diag([10 20 30]);



%% Manipulate data

csv1 = readtable('opendlv.proxy.GroundSpeedReading-112.csv');
csv2 = readtable('opendlv.proxy.AccelerationReading-112.csv');
csv3 = readtable('opendlv.proxy.AngularVelocityReading-112.csv');

vxs = csv1.groundSpeed;
vy_dots = -1*csv2.accelerationY;
phi_dots = csv3.angularVelocityZ;

vxs_temp = upsample(vxs,8);
for i = 1 : length(vxs_temp)
    vxs_temp(i) = vxs(ceil(i/8));
end
vxs = vxs_temp;
vxs = vxs(1:length(phi_dots));


%% Define Equation

X(:,1) = [0.00;0.00;0.000];
% H = [0 1 -lf; 0 0 1];
%% EKF algorithm


for i = 2:length(vxs)
    if i == 7902
        disp('here')
        
    end
    vx = vxs(i-1);
    vxc = vxs(i);
    if vx >= vx_threshold 
       F = [ 
             0           , 0                        , ((lr+lf)/vx);
             (h*cf/m)    , (1-(h*(cf+cr))/(m*vx))   , (h*(cr*lr-cf*lf)/(m*vx) - h*vx);
             (h*lf*cf/Iz), (h*(lr*cr-lf*cf)/(Iz*vx)),(1-h*(lf*lf*cf+lr*lr*cr)/(Iz*vx));    
            ];
       X_predict = F*X(:,i-1);
       P = F*P*F' + Q;
       if vxc >= vx_threshold
            H = [
                    (cf/m) , (-(cf+cf)/(m*vxc)), ((cr*lr-cf*lf)/(m*vxc) - vxc);
                    0      , 0                 , 1;
                   ];
            G = P*H'/(H*P*H' + R);
            Z = [vy_dots(i); phi_dots(i)];
            Y = Z - H*X_predict;
            X_update = X_predict + G*Y;
            P = (eye(3) - G*H)*P;
            X(:,i) = X_update;
       else
           disp('vxc is belowwww threshold')
           X(:,i) = X(:,i-1);
       end
        
    else
        disp('Vx is below threshold')
        X(:,i) = X(:,i-1);
    end
    
      
end
start_plot_at = 24200;
% start_plot_at = 13320;
range_to_plot = [start_plot_at:start_plot_at+4000];
sample_time = [0:0.005:20];
figure(1)
plot(sample_time,X(1,range_to_plot))
hold on
plot(sample_time,X(3,range_to_plot))
hold off
legend('Theta estimation','Yawrate estimation')
title('The plot of Theta and Yawrate estimated using EKF')
figure(2)
plot(sample_time,phi_dots(range_to_plot))
hold on
plot(sample_time,X(3,range_to_plot))
hold off
legend('Yawrate sensors','Yawrate estimation')
title('The plot of Yawrate from observation and estimation by EKF')
G
figure(3)
yyaxis left
plot(sample_time,X(1,range_to_plot),'-k','LineWidth',2.5)
ylim([-1 1.5])
hold on
plot(sample_time,phi_dots(range_to_plot),'-')
title('Theta estimated by EKF vs. Yawrate and Lateral acceleration measured observed during 20 seconds')

yyaxis right
ylim([-20 15])
plot(sample_time,vy_dots(range_to_plot),'--')
% legend('Theta estimation 20s','Yawrate sensors')
legend('Theta estimation 20s','Yawrate sensors','Lateral acceleration sensors')
hold off
yyaxis left
figure(4)
plot(X(1,:),'-k')

hold on
plot(X(2,:))
plot(X(3,:))
legend('Theta estimation','Yawrate estimation','Lateral speed estimation')
hold off

figure(5)
yyaxis left
plot(phi_dots,'-')
hold on
yyaxis right
plot(vy_dots,'-')
hold off
figure(6)
theta_to_plot = X(1,8500:29500);
plot(theta_to_plot,'-')
legend('Theta estimation')
title('Theta estimated from EKF')
ylim([-0.8 0.8])
xlim([1 20900])
