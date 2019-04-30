%Author: Kanadech Jirapongtanavech, Z5176970

%Program: Solution for AAS, T1 2019, Project4.part B....

function []=Project04_B()
% sim car in random movement

clc; clear all; close all; dbstop if error;
set(0,'defaultaxesfontname','times new roman');

field.range=50;
time.dt=1; time.T=500;
%% Initialization
[fig]=FigureNew(field);

carPosition = [rand*field.range,rand*field.range,wrapToPi(rand*2*pi)];
[car] = CarNew(carPosition,'r');
[car] = CarNow(car,time,0,0);
[car] = CarShow(car,0);

targetHeading = wrapToPi(rand*2*pi);

targetPosition = [rand*field.range,rand*field.range,targetHeading];
[target] = CarNew(targetPosition,'b');
[target] = CarNow(target,time,0,0);
[target] = CarShow(target,0);

d = 10;

virtualPosition = [targetPosition(1)-d*cos(targetHeading),targetPosition(2)-d*sin(targetHeading),targetHeading];
[virtual] = CarNew(virtualPosition,'k');
[virtual] = CarNow(virtual,time,0,0);
[virtual] = CarShow(virtual,0);

fis_vel=readfis('MTRN4010_vel_2.fis');
fis_ang=readfis('MTRN4010_ang.fis');
%fis_vel.DefuzzificationMethod = "centroid";
%fis_ang.DefuzzificationMethod = "mom";
%% Calculate the velocity at which the virtual target moves to real target

virtualVel = calculateDistance(virtual,target)/time.T;

%%
for t=0:time.dt:time.T
    vel = evalfis(calculateDistance(car,virtual),fis_vel);
    omega = evalfis(calculateAngle(car,virtual),fis_ang);
    
    [virtual] = CarNow(virtual,time,virtualVel,0);
    [virtual] = CarShow(virtual,t);
    
    [car]=CarNow(car,time,vel,omega);
    [car]=CarShow(car,t);
end

end

function [fig]=FigureNew(field)
fig=figure('units','normalized','position',[0.1 0.2 0.5 0.5]);
axis([-1 1 -1 1]*field.range); hold on; grid on; axis equal;
xlabel('x-direction'); ylabel('y-direction');
end

function [car]=CarNew(position,color)
car.x=position(1); car.y=position(2); car.q=position(3);
car.trace=[car.x; car.y; car.q];
car.shape=[ 2 0; 1 1; -1 1; -1 -1; 1 -1; 2 0]';
car.hdL.shape=plot(car.shape(1,:),car.shape(2,:),'color',color,'linewidth',2);
car.hdL.trace=plot(car.trace(1,:),car.trace(2,:),'color',[0 0.66 0]);
end

function [car]=CarNow(car,time,v,w)
car.x=car.x+time.dt*v*cos(car.q);
car.y=car.y+time.dt*v*sin(car.q);
car.q=car.q+time.dt*w;
car.trace(:,end+1)=[car.x; car.y; car.q];
end

function [car]=CarShow(car,t)
ax=axis;
Rz=[  cos(car.q) -sin(car.q); 
      sin(car.q)  cos(car.q)];
shape=Rz*car.shape+repmat([car.x;car.y],1,6);
set(car.hdL.shape,'xdata',shape(1,:),'ydata',shape(2,:)); 
set(car.hdL.trace,'xdata',car.trace(1,:),'ydata',car.trace(2,:));
axis(ax); title(sprintf('Time %d',t)); pause(0.001);
end

function distance = calculateDistance(src,dest)
    distance = sqrt((src.x-dest.x)^2+(src.y-dest.y)^2);
end

function omega = calculateAngle(src,dest)
    omega=atan2((dest.y-src.y),(dest.x-src.x))-src.q;
    omega=wrapToPi(omega);
end