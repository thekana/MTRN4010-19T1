A = 100; %rad/s^2
B = 2; %1/s
C = 1;
dt = 1/1000; %fix time step
xk = zeros(2,7000);
xk(:,1) = [pi/3; 0];
time = dt:dt:7;
u = (time>=1 & time<=2)*10 + (time>=3 & time<=4)*20; % better computation?
% u = zeros(1,7000);
% u(1000:2000) = 10;
% u(3000:4000) = 20;
for k = 1:size(xk,2)-1
    uk = u(k);
    xk(:,k+1) = xk(:,k) + dt*[xk(2,k); -A*sin(xk(1,k)) - B*xk(2,k) + uk];
end 
figure
plot(time,xk(1,:))
figure
plot(time,xk(2,:))