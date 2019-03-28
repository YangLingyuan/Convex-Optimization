close all;
clc; clear all;
%% Parameters of mobiler robot from the paper
vehicle.length = 0.526;  % Vehicle length
vehicle.width = 0.526;   % Vehicle width
wheel.radius = 0.0854; %
% wheel.offcentre = 0.045; %
wheel.offcentre = 0; %
a = vehicle.width/2;
b = vehicle.length/2;
rw = wheel.radius;
hx = a;
hy = b;
d = wheel.offcentre;
r = wheel.radius;

hx1 = b;
hx2 = b; 
hx3 = -b;
hx4 = -b;

hy1 = -a;
hy2 = a;
hy3 = a;
hy4 = -a;

hxi =[hx1, hx2, hx3, hx4];
hyi =[hy1, hy2, hy3, hy4];
h_i=[hxi(1),hyi(1);hxi(2),hyi(2);hxi(3),hyi(3);hxi(4),hyi(4)];
%%
R_max= 20 ;
Limit_Betadot = 2;%(rad/s)
Limit_Betadotdot = 25;%(rad/s^2)
delta1 = 1e-9;
delta2 = 1e-2;
%% Trajectory generation
%% S function defined in global coordiante system
xf=10; %m
tf=15; %s
deltaT=0.01;%s
t=0:deltaT:tf;
alpha=xf/tf;
origin_x_ICR = 1*alpha*t; %m
a=2;
c=5;
% y =(1./(1+exp(-a*(alpha*t-c))));
origin_y_ICR = 1./(1+exp(-a*(origin_x_ICR-c)));
theta = 1./(1+exp(-a*(origin_x_ICR-c)));
origin_x_ICR = origin_x_ICR;
origin_y_ICR = origin_y_ICR*0;
%% New trajectory
path_x = 1*origin_x_ICR;
path_y = 1*origin_y_ICR;

init_theta = atan((path_y(2)-path_y(1))/(path_x(2)-path_x(1)));
initPose = [path_x(1); path_y(1);init_theta];    % Initial pose (x y theta)
n = numel(path_x); % Number of grid points
theta = 1*180*theta*pi/180;

pose =[path_x; path_y; theta;];

h = mean(diff(t));
xdot = [0;diff(origin_x_ICR') / h]';
ydot = [0;diff(origin_y_ICR') / h]';
thetadot = [0;diff(theta') / h]';

velW = [xdot; ydot; thetadot]; % Velocity in global cooridante system

for k =1:n
    velB(:,k) = worldToBody(velW(:,k),pose(:,k));  % Velocity in Robot base cooridante system
end

velB = velW;
xdotB = velB(1,:);
ydotB = velB(2,:);
thetadotB = velB(3,:);
%% Acceleration in Base coordiante
xdotdotB = [0;diff(xdotB') / h]';
ydotdotB = [0;diff(ydotB') / h]';
thetadotdotB = [0;diff(thetadotB') / h]';
%
Zetadoti = [xdotB', ydotB', thetadotB'];
Zetadotdoti = [xdotdotB', ydotdotB', thetadotdotB'];
etadot =[xdotB', ydotB', thetadotB', xdotdotB', ydotdotB', thetadotdotB'];


deltaTime = t;
%% Pose store in glboal coordinate system
pose = zeros(3,n+1);
pose(:,1) = initPose;

% delta2c =[0.4444 0.4444 0.4444 0.4444];
delta2c =[0 0 0 0];
delta1 = 1e-9;


ICR_X_ref = R_max * tanh(-Zetadoti(:,2)./(Zetadoti(:,3) + delta1*sign(1-any(Zetadoti(:,3),2)+Zetadoti(:,3)))./R_max);
ICR_Y_ref = R_max * tanh(Zetadoti(:,1)./(Zetadoti(:,3) + delta1*sign(1-any(Zetadoti(:,3),2)+Zetadoti(:,3)))./R_max);
ICR_X_next = zeros(1, length(ICR_X_ref));
ICR_Y_next = zeros(1, length(ICR_Y_ref));
Beta_real=zeros(n,4);

P = [1,0;0,1];
q = [-2*ICR_X_ref(k),-2*ICR_Y_ref(k)];
r = ICR_X_ref(k)^2+ICR_Y_ref(k)^2;

Beta=zeros(n,4);
Beta_max=zeros(n,4);
Beta_min=zeros(n,4);

A=zeros(8,2);
timming=zeros(1,n);
ERR_refNext=zeros(1,n);
for k = 1:n
    if k ==1
        for i = 1:4
            yhatB = ydotB(k+1) + hxi(i) * thetadotB(k+1);
            xhatB = xdotB(k+1) - hyi(i) * thetadotB(k+1);
            Beta(k,i) = atan(yhatB/(xhatB+delta1*sign(xhatB)));
            diffBetax = -yhatB/(xhatB^2+yhatB^2+delta2c(:,i));
            diffBetay = xhatB/(xhatB^2+yhatB^2+delta2c(:,i));
            diffBetat = (hxi(i)*xdotB(k+1)+hyi(i)*ydotB(k+1))/(xhatB^2+yhatB^2+delta2c(:,i));
            F1i(i,:) = [diffBetax diffBetay diffBetat];
            F2i(i,:) = [cos(Beta(k,i)), sin(Beta(k, i)), d-hyi(i)*cos(Beta(k, i))+hxi(i)*cos(Beta(k, i))];
            feasibility=0;
            if(xhatB<0)
                if(yhatB<0)
                    Beta_real(k,i)=Beta(k,i)-pi;
                else
                    Beta_real(k,i)=Beta(k,i)+pi;
                end
            else
                Beta_real(k,i)=Beta(k,i);
            end
            
        end
    else
        for i =1:4
            yhatB = ydotB(k) + hxi(i)*thetadotB(k);
            xhatB = xdotB(k) - hyi(i)*thetadotB(k);
            
            Beta(k, i) = atan(yhatB/(xhatB+delta1*sign(xhatB)));
            diffBetax = -yhatB/(xhatB^2+yhatB^2+delta2c(:,i));
            diffBetay = xhatB/(xhatB^2+yhatB^2+delta2c(:,i));
            diffBetat = (hxi(i)*xdotB(k)+hyi(i)*ydotB(k))/(xhatB^2+yhatB^2+delta2c(:,i));
            F1i(i,:) = [diffBetax diffBetay diffBetat];
            F2i(i,:) = [cos(Beta(k, i)), sin(Beta(k, i)), d-hyi(i)*cos(Beta(k, i))+hxi(i)*sin(Beta(k, i))];
            
            if(xhatB<0)
                if(yhatB<0)
                    Beta_real(k,i)=Beta(k,i)-pi;
                else
                    Beta_real(k,i)=Beta(k,i)+pi;
                end
            else
                Beta_real(k,i)=Beta(k,i);
            end
            
        end
    end
    Betast(k,:) = Beta(k,:);
    Betaidot(k,:) = F1i * (Zetadotdoti(k,:))';
    Phi_dot(k,:) = (1/rw) * F2i * (Zetadoti(k,:))' + (d/rw) *F1i * (Zetadotdoti(k,:))';
    
    Mi = [zeros(4, 3), F1i;(1/rw)*F2i, (d/rw)*F1i];
    Abdot = (Mi*etadot(k,:)')';
    Abdotn(k,:) = (Mi*etadot(k,:)')';
    
    [ICR_next, ERR_refNext(k), timming(k), solved] = Optimization(Beta_real(k,:), Betaidot(k,:), [ICR_X_ref(k),ICR_Y_ref(k)], deltaT, h_i);
    if(solved)
        ICR_X_next(k)=ICR_next(1);
        ICR_Y_next(k)=ICR_next(2);
        ERR_refNext(k)=ICR_Y_next(k)-ICR_Y_ref(k);
    else
        ICR_X_next(k)=ICR_X_ref;
        ICR_Y_next(k)=ICR_Y_ref;
        ERR_refNext(k)=1000;
    end
    
    %% Forward kinematic model
    lamda =0.00001; % Damping factor
    F2d = pinv(F2i'*F2i+lamda^2*eye(3))*F2i';
    Mf = [-d*F2d, r*F2d];
    Zetadotnewi(k,:) = (Mf* Abdot')';
    velB(:,k) = (Mf* Abdot');
    velW(:,k) = bodyToWorld(velB(:,k),pose(:,k))';
end


Abdot = Abdotn;
allSteerAngles = Betast.*(180/pi);
allSteerRates = Betaidot.*(180/pi);
allWheelSpeeds = Phi_dot.*(180/pi);

%% Time vector separately for node and grid points
n = length(deltaTime);
timeNode = zeros(1,n+1);
timeGrid = zeros(1,n);
timeNode(1) = 0;
for i = 2:n+1
    timeNode(i) = timeNode(i-1) + deltaTime(i-1);
end
timeGrid(1) = deltaTime(1)/2;
for i = 2:n
    timeGrid(i) = deltaTime(i-1)/2 + deltaTime(i)/2;
end
%%
velWT=velW;
for k=1:n
    pose(1,k+1) = pose(1,k) + velWT(1,k)*deltaT;
    pose(2,k+1) = pose(2,k) + velWT(2,k)*deltaT;
    pose(3,k+1) = pose(3,k) + velWT(3,k)*deltaT;
end
%%
Verror = Zetadoti-Zetadotnewi;
r_path = sqrt(path_x.^2+path_y.^2);
p_path = sqrt( (pose(1,2:n+1)).^2+(pose(2,2:n+1)).^2);
p_error = r_path - p_path;

% %%
% figure(1) % trajectory
% hold on
% plot(pose(1,1),pose(2,1),'ro', ...
%      pose(1,end),pose(2,end),'go', ...
%      pose(1,:),pose(2,:),'r-');
% axis equal
% title('Vehicle Trajectory');
% xlabel('X [m]')
% ylabel('Y [m]')
% hold on
% plot(path_x,path_y)
% legend('Start','End','Trajectory', 'Original Path')
%
% %%
% figure(2) % trajectory
% plot(timeGrid,p_error)
% xlabel('time [s]')
% ylabel('Error [m]')
% title('Trajectory path error');
% legend('Path Error [m]')

figure(2) % trajectory
plot(timeGrid,ERR_refNext);
hold on;
plot(timeGrid,ICR_Y_next);

hold on;
plot(timeGrid,ICR_Y_ref);

xlabel('time [s]')
ylabel('Error [m]')
title('ICR Error ref-curr');
legend('ICR Error [m]','ICR-Y-next','ICR-Y-ref')

figure(3)
plot(timeGrid,timming)
xlabel('time [s]')
ylabel('optimization timming [s]')
title('optimization timming per iteration');
legend('opt timming [m]')
%
% figure(3) % steering angles
% plot(timeGrid,allSteerAngles(:,1),...
%      timeGrid,allSteerAngles(:,2),...
%      timeGrid,allSteerAngles(:,3),...
%      timeGrid,allSteerAngles(:,4));
% legend('Wheel 1','Wheel 2','Wheel 3','Wheel 4')
% title('Wheel Steering Angles');
% xlabel('Time [s]');
% ylabel('Wheel steering angle [Deg]');
%
%
% figure(4) % steering angles
% plot(timeGrid,allSteerRates(:,1),...
%      timeGrid,allSteerRates(:,2),...
%      timeGrid,allSteerRates(:,3),...
%      timeGrid,allSteerRates(:,4));
% legend('Wheel 1','Wheel 2','Wheel 3','Wheel 4')
% title('Wheel Steering Rates');
% xlabel('Time [s]');
% ylabel('Wheel steering angle [Deg/s]');
%
%
% figure(5) % wheel angular speeds
% plot(timeGrid,allWheelSpeeds(:,1),...
%      timeGrid,allWheelSpeeds(:,2),...
%      timeGrid,allWheelSpeeds(:,3),...
%      timeGrid,allWheelSpeeds(:,4));
% legend('Wheel 1','Wheel 2','Wheel 3','Wheel 4')
% title('Wheel Angular Speeds');
% xlabel('Time [s]');
% ylabel('Wheel angular speed [Deg/s]');
%
%
% figure(6)
% subplot(3,1,1)
% plot([Zetadoti(:,1),Zetadotnewi(:,1)]);
% legend('Vxapplied', 'Vxforward')
% xlabel('time (sec)')
% ylabel('Platfrom Speed Vx(rad/s))')
% title('Platform Vx');
% subplot(3,1,2)
% plot([Zetadoti(:,2),Zetadotnewi(:,2)]);
% legend('Vyapplied', 'Vyforward')
% xlabel('time (sec)')
% ylabel('Platfrom Speed Vy(rad/s))')
% title('Platform Vy');
% subplot(3,1,3)
% plot([Zetadoti(:,3)]);
% hold on;
% plot([Zetadotnewi(:,3)], '--');
% legend('Vzapplied', 'Vzforward')
% xlabel('time (sec)')
% ylabel('Rotational speed \theta(rad/s))')
% title('Platfrom rotationanl speed');
%
% figure(7)
% plot([Verror(:,1),Verror(:,2), Verror(:,3)]);
% legend('X_{Error}', 'Y_{Error}', 'Z_{Error}')
% title('Velocity error beteen applied and calculated')
%
% figure(8)
% plot([Zetadoti(:,2),Zetadotnewi(:,2)]);
% legend('Orignal', 'Analytical')
%
% figure(9)
% plot(t, velWT(1,:))
% hold on;
% plot(t, velWT(2,:))
% title('Vel')

%% Cube rotation and translation
% New plot
origin_x_ICR = [-0.1 -0.1 0.1 0.1];
origin_y_ICR = [-0.1 0.1 0.1 -0.1];
origin_x_body = [hx4 hx3 hx2 hx1]; %// initial coordinates of vertices
origin_y_body = [hy4 hy3 hy2 hy1];
origin_x = [origin_x_ICR(1) origin_x_body(1); origin_x_ICR(2) origin_x_body(2); origin_x_ICR(3) origin_x_body(3); origin_x_ICR(4) origin_x_body(4)]; %// initial coordinates of vertices
origin_y = [origin_y_ICR(1) origin_y_body(1); origin_y_ICR(2) origin_y_body(2); origin_y_ICR(3) origin_y_body(3); origin_y_ICR(4) origin_y_body(4)];
c = [0;1];
destination_x = origin_x_body + 0.25; %// final coordinates of vertices
destination_y = origin_y_body + 0.25;
n_steps = 100; %// number of "frames"
t_pause = 0.1; %// seconds between frames





figure(10)

h = patch(origin_x, origin_y, c); %// create object at initial position
% axis([-1 15 -1 1]) %// adjust as needed, to cover the desired area
xlim([-1 10])
ylim([-10 1])
%axis equal %// same scale in both axes
axis manual %// prevent axes from auto-scaling

hb = [origin_x_body; origin_y_body];
cb = [origin_x_ICR; origin_y_ICR];
count=1;
theta = pose(3,:);

thetanew = [theta, theta+pi];

for i=1:25:length(theta)
    R=[cos(theta(i)) sin(theta(i));-sin(theta(i)) cos(theta(i))];
    hi = R*hb;
    ci = R*cb;
    x1 = pose(1,i)+hi(1,:);
    y1 = pose(2,i)+hi(2,:);
    x2 = pose(1,i)+ICR_X_ref(i)+ci(1,:);
    y2 = pose(2,i)+ICR_Y_ref(i)+ci(2,:);
    
    patch_pose_x = [x1',x2']; %// initial coordinates of vertices;
    patch_pose_y = [y1',y2'];
    h.XData = patch_pose_x;
    h.YData = patch_pose_y;
    %set(h, 'Vertices', [x1(:) y1(:)]) %// change object's position
    
    pause(t_pause) %// a pause is needed to make movement slower
    drawnow %// probably not needed after pause. Just in case
end

%% Animated plots

stangles = [allSteerAngles; allSteerAngles];
SteeRates = [allSteerRates; allSteerRates];
WheelSpeeds = [allWheelSpeeds; allWheelSpeeds];


for i=1:length(stangles)-1
    stanglesnew(i,:)=stangles(i,:);
    SteeRatesnew(i,:)=SteeRates(i,:);
    WheelSpeedsnew(i,:)=WheelSpeeds(i,:);
end
%
% figure(12)
%     h = animatedline('Color','b');
%     h2 = animatedline('Color','r');
%     h3 = animatedline('Color','g');
%     h4 = animatedline('Color','k');
%
%     xlabel('time (sec)')
%     ylabel('Steering angles (Degrees))')
%     timeGrid = t;
%     allSteerAngles = stanglesnew;
%
% numpoints = length(timeGrid);
%
%     for k = 1:25:numpoints
%       addpoints(h,timeGrid(k),allSteerAngles(k,1))
%       addpoints(h2,timeGrid(k),allSteerAngles(k,2))
%       addpoints(h3,timeGrid(k),allSteerAngles(k,3))
%       addpoints(h4,timeGrid(k),allSteerAngles(k,4))
%       legend('Steering angle1', 'Steering angle2', 'Steering angle3', 'Steering angle4')
%       drawnow update;
%       frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1;
%       imwrite(imind,cm,'SteerAngle.gif','gif','Loopcount',inf,'DelayTime',0.01);
%     else
%       imwrite(imind,cm,'SteerAngle.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
%     end
%
%     figure(13)
%     h = animatedline('Color','b');
%     h2 = animatedline('Color','r');
%     h3 = animatedline('Color','g');
%     h4 = animatedline('Color','k');
%
%     xlabel('time (sec)')
%     ylabel('Steering Rates (Degrees/sec))')
%     timeGrid = t;
%     allSteerRates = SteeRatesnew;
%
% numpoints = length(timeGrid);
%
%     for k = 1:25:numpoints
%       addpoints(h,timeGrid(k),allSteerRates(k,1))
%       addpoints(h2,timeGrid(k),allSteerRates(k,2))
%       addpoints(h3,timeGrid(k),allSteerRates(k,3))
%       addpoints(h4,timeGrid(k),allSteerRates(k,4))
%       legend('Steering Rate1', 'Steering Rate2', 'Steering Rate3', 'Steering Rate4')
%       drawnow update;
%       frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1;
%       imwrite(imind,cm,'SteeRate.gif','gif','Loopcount',inf,'DelayTime',0.01);
%     else
%       imwrite(imind,cm,'SteeRate.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
%     end
%     %%
%     figure(14)
%     h = animatedline('Color','b');
%     h2 = animatedline('Color','r');
%     h3 = animatedline('Color','g');
%     h4 = animatedline('Color','k');
%
%     xlabel('time (sec)')
%     ylabel('Wheel Speed (Degrees/sec))')
%     timeGrid = t;
%     allWheelSpeeds = WheelSpeedsnew;
%
% numpoints = length(timeGrid);
%
%     for k = 1:25:numpoints
%       addpoints(h,timeGrid(k),allWheelSpeeds(k,1))
%       addpoints(h2,timeGrid(k),allWheelSpeeds(k,2))
%       addpoints(h3,timeGrid(k),allWheelSpeeds(k,3))
%       addpoints(h4,timeGrid(k),allWheelSpeeds(k,4))
%       legend('Wheel Speed1', 'Wheel Speed2', 'Wheel Speed3', 'Wheel Speed4')
%       drawnow update;
%       frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1;
%       imwrite(imind,cm,'WheelSpeed.gif','gif','Loopcount',inf,'DelayTime',0.01);
%     else
%       imwrite(imind,cm,'WheelSpeed.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
%     end
%
%     %% Combined animation plot with synchronization
% figure(15)
%     subplot(2,2,1:2)
%     % plot(T, Y, 'LineWidth', 2)
%     h = fill(origin_x_body, origin_y_body, 'r'); %// create object at initial position
%     xlabel('Distance in X-axis (meters)')
%     ylabel('Distance in Y-axis (meters)')
%     xlim([-1 15])
%     ylim([-1 1])
%     axis equal
%     subplot(2,2,3)
%     h1 = animatedline('Color','b');
%     h2 = animatedline('Color','r');
%     h3 = animatedline('Color','g');
%     h4 = animatedline('Color','k');
%     legend('ST1', 'ST2', 'ST3', 'ST4')
%
%     xlabel('time (sec)')
%     ylabel('Steering angles (Degrees))')
%
%     subplot(2,2,4)
%     h5 = animatedline('Color','b');
%     h6 = animatedline('Color','r');
%     h7 = animatedline('Color','g');
%     h8 = animatedline('Color','k');
%     legend('WS1', 'WS2', 'WS3', 'WS4')
%
%     xlabel('time (sec)')
%     ylabel('Wheel Speed (Degrees/sec))')
%
%     timeGrid = t;
%     allSteerAngles = stanglesnew;
%
% for i=1:25:length(theta)
%     R=[cos(theta(i)) sin(theta(i));-sin(theta(i)) cos(theta(i))];
%     hi = R*hb;
% %     x = path_x(i)+hi(1,:);
% %     y = path_y(i)+hi(2,:);
%
%     origin_x_ICR = pose(1,i)+hi(1,:);
%     origin_y_ICR = pose(2,i)+hi(2,:);
%      set(h, 'Vertices', [origin_x_ICR(:) origin_y_ICR(:)]) %// change object's position
%      addpoints(h1,timeGrid(i),allSteerAngles(i,1))
%      addpoints(h2,timeGrid(i),allSteerAngles(i,2))
%      addpoints(h3,timeGrid(i),allSteerAngles(i,3))
%      addpoints(h4,timeGrid(i),allSteerAngles(i,4))
% %      legend('ST1', 'ST2', 'ST3', 'ST4')
%      drawnow update;
%      addpoints(h5,timeGrid(i),allWheelSpeeds(i,1))
%       addpoints(h6,timeGrid(i),allWheelSpeeds(i,2))
%       addpoints(h7,timeGrid(i),allWheelSpeeds(i,3))
%       addpoints(h8,timeGrid(i),allWheelSpeeds(i,4))
% %       legend('WS1', 'WS2', 'WS3', 'WS4')
%       drawnow update;
%
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 1;
%       imwrite(imind,cm,'platformcom.gif','gif','Loopcount',inf,'DelayTime',0.01);
%     else
%       imwrite(imind,cm,'platformcom.gif','gif','WriteMode','append','DelayTime',0.01);
%     end
% end

%% ICR Based control of mobile robot
Beta1 = 70.36;
Beta2 = 109.64;
Beta3 = 101.76;
Beta4 = 78.24;

B12 = Beta1-Beta2;
B13 = Beta1-Beta3;
B14 = Beta1-Beta4;
B23 = Beta2-Beta3;
B24 = Beta2-Beta4;
B34 = Beta3-Beta4;

Bij = abs([B12, B13, B14, B23, B24, B34])*pi/180;

% for i=1:length(Bij)
%     if Bij(i)<=(pi/2)
%         elseif Bij(i)>(pi/2)
%             Bij(i) =Bij(i)-pi;
%     end
% end




Betai = Beta1*pi/180;
Betaj = Beta2*pi/180;
Beta_th = 80*pi/180;
R = 10;

if Betai < Beta_th
    Ycur = R *cos(Betai);
    Xcur = R *sin(Betai);
else
    Ycur = (hxi(1)-hxi(2)+hyi(1)*tan(Betai)-hyi(2)*tan(Betaj))/(tan(Betai)-tan(Betaj));
    Xcur = hxi(1)-(Ycur-hyi(1))*tan(Betai);
end



