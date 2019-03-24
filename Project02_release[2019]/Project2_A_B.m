load('DataForProject02/IMU_dataC.mat');
load('DataForProject02/Speed_dataC.mat');
load('DataForProject02/Laser__2C.mat');

% remove bias
time = double(IMU.times-IMU.times(1))/10000;
Laser_time = double(dataL.times-dataL.times(1))/10000;

yaw = -IMU.DATAf(6,:); % negate w(K)

horizon0 = length(find(time<20)); %stationary for 20 secs
Bias = mean(yaw(1:horizon0)); %Bias

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

thetaKL = zeros(1,length(Laser_time));
xKL = zeros(1,length(Laser_time));
yKL = zeros(1,length(Laser_time));
thetaKL(1) = pi/2;
xKL(1) = 0;
yKL(1) = 0;

j = 2;
for i = 2:L
    dt = time(i)-time(i-1); % find dt
    thetaK(i) = thetaK(i-1) + dt * yawC(i-1); % euler approx
    xK(i) = xK(i-1) + dt * Vel.speeds(i-1)*cos(thetaK(i));
    yK(i) = yK(i-1) + dt * Vel.speeds(i-1)*sin(thetaK(i));
    if (j <= length(Laser_time) && Laser_time(j) - time(i-1)< dt)
        dtL = Laser_time(j) - time(i-1);
        thetaKL(j) = thetaK(i-1) + dtL * yawC(i-1); % euler approx
        xKL(j) = xK(i-1) + dtL * Vel.speeds(i-1)*cos(thetaK(i));
        yKL(j) = yK(i-1) + dtL * Vel.speeds(i-1)*sin(thetaK(i));
        j = j + 1;
    end
end

% convert from radian to degree
thetaK = thetaK * 180/pi;
thetaKL = thetaKL * 180/pi;
figure();
plot(time,thetaK);
figure();
plot(xK,yK);
figure();
plot(Laser_time,thetaKL);
figure();
plot(xKL,yKL);