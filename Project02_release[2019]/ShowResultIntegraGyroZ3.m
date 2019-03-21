
% axample about inegrating Wz ( we assume Wx=Wy=0; for our problem, in which
% the platform is assumed to have pitch=roll=0; in Project2

% MTRN4010, 2019.T1



% load provided IMU data
load('DataForProject02\IMU_dataC.mat');
k=180/pi;

wz=-IMU.DATAf(6,:);             %  GyroZ, in rad/sec
% according to my coordinate frame convention, I change the sign of Wz; you
% may or may not need to do it.

t0 = IMU.times(1);            % time0, first sample time;

% get average of initial samples, because the machine was not moving during
% that period of time; we take advantage of that fact, so we estimate the bias using this OFF_LINE approach  / trick.
horizon0  = 1500; % this is about  1500*5ms = 7.5 secs.
wzBias = mean(wz(1:horizon0));  
wzC = wz-wzBias;        % remove bias, 

L=length(wzC);

% create some buffers, for recording variables' contents.
aaC=zeros(1,L);
aa =zeros(1,L);

%a0 = 90*pi/180;
a0 = 0;     % initial angle. % you may use a different one.
aC=a0; a=a0;

ts= double(IMU.times-t0)*0.0001;    % time in seconds;

% I implement this in a loop, presuming you may add extra code, for
% implementing full kinematic model.
%otherwise, I would preffer to use Matlab muilt-in vector operations, which are more efficient. 
% For the full kinematic model, however, we need to use the loop.

aa(1)=a;
aaC(1)=aC;

for i=2:L;
    dt=ts(i)-ts(i-1);
    a=a+dt*wz(i-1);                        % integrating raw data 
    aC=aC+dt*wzC(i-1);                 % integrating improved  data (bias mitigated)
    
    % record in buffers, for plotting,etc.
    aa(i)=a;    
    aaC(i) = aC;
    
end;

% done.


%----------------------------------------------------------
% just plotting results.
figure(1) ; clf() ; 

subplot(211); plot(k*wz);  grid on; ax=axis(); title(' yaw rate (raw)'); ylabel(' deg/sec');
subplot(212); plot(k*wzC);grid on; axis(ax); title(' yaw rate (bias removed)'); ylabel(' deg/sec');

figure(2) ; clf();

plot(ts,k*aaC,'r');hold on;
plot(ts,k*aa,'b'); 
xlabel('time (in seconds)');

legend({'corrected','biased'});
grid on;
title(' yaw rate integrated (assuming yaw0=0)');
ylabel('yaw (degrees)');
%----------------------------------------------------------




