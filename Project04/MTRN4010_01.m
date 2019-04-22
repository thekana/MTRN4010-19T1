function []=MTRN4010_01()
% sim car in random movement

clc; clear all; close all; dbstop if error;
set(0,'defaultaxesfontname','times new roman');

field.range=50;
time.dt=1; time.T=500;

[fig]=FigureNew(field);
[car]=CarNew();

for t=0:time.dt:time.T,
  [car]=CarNow(car,time,rand,randn*0.33);
  [car]=CarShow(car,t);
end;


function [fig]=FigureNew(field)
fig=figure('units','normalized','position',[0.1 0.2 0.5 0.5]);
axis([-1 1 -1 1]*field.range); hold on; grid on; axis equal;
xlabel('x-direction'); ylabel('y-direction');

function [car]=CarNew()
car.x=0; car.y=0; car.q=0;
car.trace=[car.x; car.y; car.q];
car.shape=[ 2 0; 1 1; -1 1; -1 -1; 1 -1; 2 0]';
car.hdL.shape=plot(car.shape(1,:),car.shape(2,:),'color','b','linewidth',2);
car.hdL.trace=plot(car.trace(1,:),car.trace(2,:),'color',[0 0.66 0]);

function [car]=CarNow(car,time,v,w)
car.x=car.x+time.dt*v*cos(car.q);
car.y=car.y+time.dt*v*sin(car.q);
car.q=car.q+time.dt*w;
car.trace(:,end+1)=[car.x; car.y; car.q];

function [car]=CarShow(car,t)
ax=axis;
Rz=[  cos(car.q) -sin(car.q); 
      sin(car.q)  cos(car.q)];
shape=Rz*car.shape+repmat([car.x;car.y],1,6);
set(car.hdL.shape,'xdata',shape(1,:),'ydata',shape(2,:)); 
set(car.hdL.trace,'xdata',car.trace(1,:),'ydata',car.trace(2,:));
axis(ax); title(sprintf('Time %d',t)); pause(0.001);

