load('DataForProject02/IMU_dataC.mat');
load('DataForProject02/Speed_dataC.mat');

% remove bias
time = double(IMU.times-IMU.times(1))/10000;
%timeV = double(Vel.times-Vel.times(1))/10000;

yaw = -IMU.DATAf(6,:); % negate w(K)

horizon0 = length(find(time<20)); %stationary for 20 secs
Bias = mean(yaw(1:4000)); %Bias

yawC = yaw - Bias; % w(K) corrected
L = length(yawC);

% buffer for variables
thetaK = zeros(1,L);
xK = zeros(1,L);
yK = zeros(1,L);
% setting initial conditions
thetaK(1) = pi/2;
xK(1) = 0;
yK(1) = 0;

for i = 2:L
    dt = time(i)-time(i-1); % find dt
    thetaK(i) = thetaK(i-1) + dt * yawC(i-1); % euler approx
    xK(i) = xK(i-1) + dt * Vel.speeds(i-1)*cos(thetaK(i));
    yK(i) = yK(i-1) + dt * Vel.speeds(i-1)*sin(thetaK(i));
end

% convert from radian to degree
thetaK = thetaK * 180/pi;
figure();
plot(time,thetaK);
figure();
plot(xK,yK);
