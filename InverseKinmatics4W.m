function [initialPoseBus,  platformControlBus, kinematicControlFault] = InverseKinmatics4W(enable, setInitialPosition, ref, refDotDot, maximumSpeed, sampleTimeController)
% File initialization

persistent initialPositionIsSet;
if isempty(initialPositionIsSet)
    initialPositionIsSet = false;
end

persistent previousBeta;
if isempty(previousBeta)
    previousBeta = zeros(4,1);
end

persistent fault;
if isempty(fault)
    fault = false;
end
vehicle.length = 0.526;  % Vehicle length
vehicle.width = 0.526;   % Vehicle width
wheel.radius = 0.084; %
wheel.offcentre = 0.0; %
a = vehicle.width/2;
b = vehicle.length/2;
rw = wheel.radius;
hx = b;
hy = a;
d = wheel.offcentre;
%allWheelSpeeds=zeros(4,1);
hx1 = hx;
hx2 = hx;
hx3 = -hx;
hx4 = -hx;

hy1 = -hy;
hy2 =  hy;
hy3 =  hy;
hy4 = -hy;

hxi =[hx1, hx2, hx3, hx4];
hyi =[hy1, hy2, hy3, hy4];

%
%Vx = ref.xDot;
%Vy = ref.yDot;
%omega = ref.thetaDot;

%Ax = refDotDot(1);
%Ay = refDotDot(2);
%Aomega = refDotDot(3);

% Trajectory generation
%xdot = Vx;
%ydot = Vy;
%thetadot = omega;

%velW = [xdot; ydot; thetadot]; % Velocity in global cooridante system

% Acceleration in Base coordiante
%xdotdotW = Ax;
%ydotdotW = Ay;
%thetadotdotW = Aomega;

%AccW = [xdotdotW; ydotdotW; thetadotdotW]; % Velocity in global cooridante system

% Velocity and accelartion in Base coordiante system
%velB = velW;
%AccB = AccW;

%xdotB = velB(1);
%ydotB = velB(2);
%thetadotB = velB(3);

%xdotdotB = AccB(1);
%ydotdotB = AccB(2);
%thetadotdotB = AccB(3);

%
Zetadoti = [ref.xDot, ref.yDot, ref.thetaDot];
Zetadotdoti = [refDotDot(1), refDotDot(2), refDotDot(3)];

%Zetadoti = [xdotB, ydotB, thetadotB];
%Zetadotdoti = [xdotdotB, ydotdotB, thetadotdotB];

delta2c =[0.0 0.0 0.0 0.0];
delta2c =[0.4726 0.7527 0.4736 0.7518];

delta1 = 1e-9;

Beta = zeros(4,1);
Phi_dot = zeros(4,1);
Betaidot = zeros(4,1);

%position = zeros(4);
%position(2) = inputBus.steerMotor1Msg.actualPosition;
%position(1) = inputBus.steerMotor2Msg.actualPosition;
%position(4) = inputBus.steerMotor3Msg.actualPosition;
%position(3) = inputBus.steerMotor4Msg.actualPosition;
if ~enable && ~setInitialPosition
    initialPositionIsSet = false;
end

if ~enable && ~initialPositionIsSet
    previousBeta = zeros(4,1);
    fault = false;
end

% Pose store in glboal coordinate system
%if  ~isnan(xdotB) && ~isnan(ydotB) && ~isnan(thetadotB) && ~isnan(xdotdotB) && ~isnan(ydotdotB) && ~isnan(thetadotdotB)
if  ~isnan(Zetadoti(1)) && ~isnan(Zetadoti(2)) && ~isnan(Zetadoti(3)) && ~isnan(Zetadotdoti(1)) && ~isnan(Zetadotdoti(2)) && ~isnan(Zetadotdoti(3))
    if (ref.xDot ~= 0 || ref.yDot ~= 0 || ref.thetaDot ~= 0)
        for i = 1:4
            %Betai = atan((ydotB+hxi(i)*thetadotB)/(xdotB-hyi(i)*thetadotB));
            %F1i = [(-hxi(i) * thetadotB - ydotB) / ((hxi(i) ^ 2 + hyi(i) ^ 2) * thetadotB ^ 2 + (2 * hxi(i) * ydotB - 2 * hyi(i) * xdotB) * thetadotB + xdotB ^ 2 + ydotB ^ 2) (-hyi(i) * thetadotB + xdotB) / ((hxi(i) ^ 2 + hyi(i) ^ 2) * thetadotB ^ 2 + (2 * hxi(i) * ydotB - 2 * hyi(i) * xdotB) * thetadotB + xdotB ^ 2 + ydotB ^ 2) (hxi(i) * xdotB + hyi(i) * ydotB) / ((hxi(i) ^ 2 + hyi(i) ^ 2) * thetadotB ^ 2 + (2 * hxi(i) * ydotB - 2 * hyi(i) * xdotB) * thetadotB + xdotB ^ 2 + ydotB ^ 2);];
            %F2i = [cos(Betai), sin(Betai), d-hyi(i)*cos(Betai)+hxi(i)*cos(Betai)];
            %Betaist(i,:) = atan((ydotB+hxi(i)*thetadotB)/(xdotB-hyi(i)*thetadotB));
            yhatB = ref.yDot + hxi(i)*ref.thetaDot;
            xhatB = ref.xDot - hyi(i)*ref.thetaDot;
            
            Beta(i) = atan(yhatB/(xhatB+delta1*sign(xhatB)));
            diffBetax = -yhatB/(xhatB^2+yhatB^2+delta2c(i));
            diffBetay = xhatB/(xhatB^2+yhatB^2+delta2c(i));
            diffBetat = (hxi(i)*ref.xDot+hyi(i)*ref.yDot)/(xhatB^2+yhatB^2+delta2c(i));
            F1i = [diffBetax diffBetay diffBetat];
            F2i = [cos(Beta(i)), sin(Beta(i)), d-hyi(i)*cos(Beta(i))+hxi(i)*sin(Beta(i))];
            
            if abs(Beta(i) - previousBeta(i)) > pi + pi/2
                if previousBeta(i) > 0
                    Beta(i) = -sign(Beta(i))*Beta(i) + 2*pi;
                    Betaidot(i) = F1i * Zetadotdoti';
                    Phi_dot(i) = (1/rw) * F2i * Zetadoti' + (d/rw) * F1i * Zetadotdoti';
                else
                    Beta(i) = sign(Beta(i))*Beta(i) - 2*pi;
                    Betaidot(i) = F1i * Zetadotdoti';
                    Phi_dot(i) = (1/rw) * F2i * Zetadoti' + (d/rw) * F1i * Zetadotdoti';
                end
            elseif abs(Beta(i) - previousBeta(i)) > pi/2
                if previousBeta(i) > 0
                    Beta(i) = Beta(i) + pi;
                    Betaidot(i) = F1i * Zetadotdoti';
                    Phi_dot(i) = -  (1/rw) * F2i * Zetadoti' + (d/rw) * F1i * Zetadotdoti';
                else
                    Beta(i) = Beta(i) - pi;
                    Betaidot(i) = F1i * Zetadotdoti';
                    Phi_dot(i) = -(1/rw) * F2i * Zetadoti' + (d/rw) * F1i * Zetadotdoti';
                end
            else
                %Beta(i,1) = sign(Beta(i,1))*Beta(i,1);
                Betaidot(i) = F1i * Zetadotdoti';
                Phi_dot(i) = (1/rw) * F2i * Zetadoti' + (d/rw) * F1i * Zetadotdoti';
            end
        end
        
        if enable || (setInitialPosition && ~initialPositionIsSet)
            previousBeta = Beta;
            initialPositionIsSet = true;
        end
        
        % elseif( abs(previousBeta(1) - inputBus.steerMotor2Msg.actualPosition) < maximumSpeed*4*sampleTimeController &&...
        %         abs(previousBeta(2) - inputBus.steerMotor1Msg.actualPosition) < maximumSpeed*4*sampleTimeController &&...
        %         abs(previousBeta(3) - inputBus.steerMotor4Msg.actualPosition) < maximumSpeed*4*sampleTimeController &&...
        %         abs(previousBeta(4) - inputBus.steerMotor3Msg.actualPosition) < maximumSpeed*4*sampleTimeController)
    elseif ref.xDot == 0 && ref.yDot == 0 && ref.thetaDot == 0 && (enable || setInitialPosition || initialPositionIsSet)
        if setInitialPosition && ~initialPositionIsSet %previous position is keept at 0 0 0           
            initialPositionIsSet = true;
        end
        Beta = previousBeta;
    end
else
    fault = true;
end

if isnan(Beta(1)) || isnan(Betaidot(1)) || isnan(Phi_dot(1)) ||...
        isnan(Beta(2)) || isnan(Betaidot(2)) || isnan(Phi_dot(2)) ||...
        isnan(Beta(3)) || isnan(Betaidot(3)) || isnan(Phi_dot(3)) ||...
        isnan(Beta(4)) || isnan(Betaidot(4)) || isnan(Phi_dot(4))
    fault = true;
    initialPoseBus.SteerMotor1Position = 0;
    initialPoseBus.SteerMotor2Position = 0;
    initialPoseBus.SteerMotor3Position = 0;
    initialPoseBus.SteerMotor4Position = 0;
else
    initialPoseBus.SteerMotor1Position = Beta(2);
    initialPoseBus.SteerMotor2Position = Beta(1);
    initialPoseBus.SteerMotor3Position = Beta(4);
    initialPoseBus.SteerMotor4Position = Beta(3);
end

%if  enable && (...
%        abs(Beta(1) - inputBus.steerMotor2Msg.actualPosition) > maximumSpeed*sampleTimeController*4 ||...
%        abs(Beta(2) - inputBus.steerMotor1Msg.actualPosition) > maximumSpeed*sampleTimeController*4 ||...
%        abs(Beta(3) - inputBus.steerMotor4Msg.actualPosition) > maximumSpeed*sampleTimeController*4 ||...
%        abs(Beta(4) - inputBus.steerMotor3Msg.actualPosition) > maximumSpeed*sampleTimeController*4)
%    fault = true;
%end

if enable && ~fault    
    platformControlBus.steerMotor1ControlMsg.command = CANCommands.CC_PositionControl;
    platformControlBus.steerMotor1ControlMsg.position = Beta(2);
    platformControlBus.steerMotor1ControlMsg.speed = single(Betaidot(2));
    platformControlBus.steerMotor1ControlMsg.acceleration = single(0);
    
    platformControlBus.steerMotor2ControlMsg.command = CANCommands.CC_PositionControl;
    platformControlBus.steerMotor2ControlMsg.position = Beta(1);
    platformControlBus.steerMotor2ControlMsg.speed = single(Betaidot(1));
    platformControlBus.steerMotor2ControlMsg.acceleration = single(0);
    
    platformControlBus.steerMotor3ControlMsg.command = CANCommands.CC_PositionControl;
    platformControlBus.steerMotor3ControlMsg.position = Beta(4);
    platformControlBus.steerMotor3ControlMsg.speed = single(Betaidot(4));
    platformControlBus.steerMotor3ControlMsg.acceleration = single(0);
    
    platformControlBus.steerMotor4ControlMsg.command = CANCommands.CC_PositionControl;
    platformControlBus.steerMotor4ControlMsg.position = Beta(3);
    platformControlBus.steerMotor4ControlMsg.speed = single(Betaidot(3));
    platformControlBus.steerMotor4ControlMsg.acceleration = single(0);
    
    platformControlBus.wheelMotor1ControlMsg.command = CANCommands.CC_SpeedControl;
    platformControlBus.wheelMotor1ControlMsg.position = 0;
    platformControlBus.wheelMotor1ControlMsg.speed = single(Phi_dot(2));
    platformControlBus.wheelMotor1ControlMsg.acceleration = single(0);
    
    platformControlBus.wheelMotor2ControlMsg.command = CANCommands.CC_SpeedControl;
    platformControlBus.wheelMotor2ControlMsg.position = 0;
    platformControlBus.wheelMotor2ControlMsg.speed = single(Phi_dot(1));
    platformControlBus.wheelMotor2ControlMsg.acceleration = single(0);
    
    platformControlBus.wheelMotor3ControlMsg.command = CANCommands.CC_SpeedControl;
    platformControlBus.wheelMotor3ControlMsg.position = 0;
    platformControlBus.wheelMotor3ControlMsg.speed = single(Phi_dot(4));
    platformControlBus.wheelMotor3ControlMsg.acceleration = single(0);
    
    platformControlBus.wheelMotor4ControlMsg.command = CANCommands.CC_SpeedControl;
    platformControlBus.wheelMotor4ControlMsg.position = 0;
    platformControlBus.wheelMotor4ControlMsg.speed = single(Phi_dot(3));
    platformControlBus.wheelMotor4ControlMsg.acceleration = single(0);
    
else

    
    platformControlBus.steerMotor1ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.steerMotor1ControlMsg.position = 0;
    platformControlBus.steerMotor1ControlMsg.speed = single(0);
    platformControlBus.steerMotor1ControlMsg.acceleration = single(0);
    
    platformControlBus.steerMotor2ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.steerMotor2ControlMsg.position = 0;
    platformControlBus.steerMotor2ControlMsg.speed = single(0);
    platformControlBus.steerMotor2ControlMsg.acceleration = single(0);
    
    platformControlBus.steerMotor3ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.steerMotor3ControlMsg.position = 0;
    platformControlBus.steerMotor3ControlMsg.speed = single(0);
    platformControlBus.steerMotor3ControlMsg.acceleration = single(0);
    
    platformControlBus.steerMotor4ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.steerMotor4ControlMsg.position = 0;
    platformControlBus.steerMotor4ControlMsg.speed = single(0);
    platformControlBus.steerMotor4ControlMsg.acceleration = single(0);
    
    
    platformControlBus.wheelMotor1ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.wheelMotor1ControlMsg.position = 0;
    platformControlBus.wheelMotor1ControlMsg.speed = single(0);
    platformControlBus.wheelMotor1ControlMsg.acceleration = single(0);
    
    platformControlBus.wheelMotor2ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.wheelMotor2ControlMsg.position = 0;
    platformControlBus.wheelMotor2ControlMsg.speed = single(0);
    platformControlBus.wheelMotor2ControlMsg.acceleration = single(0);
    
    platformControlBus.wheelMotor3ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.wheelMotor3ControlMsg.position = 0;
    platformControlBus.wheelMotor3ControlMsg.speed = single(0);
    platformControlBus.wheelMotor3ControlMsg.acceleration = single(0);
    
    platformControlBus.wheelMotor4ControlMsg.command = CANCommands.CC_Stop;
    platformControlBus.wheelMotor4ControlMsg.position = 0;
    platformControlBus.wheelMotor4ControlMsg.speed = single(0);
    platformControlBus.wheelMotor4ControlMsg.acceleration = single(0);
end

kinematicControlFault = fault;
%allSteerSpeeds = Betaidot;
%allSteerAngles = Beta;
%allWheelSpeeds = Phi_dot;
