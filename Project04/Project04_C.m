%Author: Kanadech Jirapongtanavech, Z5176970

%Program: Solution for AAS, T1 2019, Project4.part C....

function []=Project04_C()
% sim car in random movement

clc; clear ; close all; dbstop if error;
set(0,'defaultaxesfontname','times new roman');
field.range=50;
time.dt=1; time.T=500;
%% Initialization

[fig] = FigureNew(field);
carPosition = [(rand-0.5)*field.range,(rand-0.5)*field.range,wrapToPi(rand*2*pi)];
[car] = CarNew(carPosition,'r');
[car] = CarNow(car,time,0,0);

targetHeading = wrapToPi(rand*2*pi);

targetPosition = [(rand-0.5)*field.range,(rand-0.5)*field.range,targetHeading];
[target] = CarNew(targetPosition,'b');
[target] = CarNow(target,time,0,0);


%% import fuzzy system
fis_vel=readfis('MTRN4010_vel_2.fis');
fis_ang=readfis('MTRN4010_ang.fis');

%% PSO parameters

PSO.DLB=0; PSO.DUB=30; %distance lower bound and upper bound
PSO.D=1; PSO.G=20; PSO.N=20;% particle dimension generations, particles
PSO.V=rand(PSO.D,PSO.N);% initial PSO particle velocity
PSO.Gbest=[]; PSO.gbest=realmax;% PSO gbest
PSO.Pbest=[]; PSO.pbest=ones(1,PSO.N)*realmax;% PSO pbest
PSO.w1=0.9; PSO.w2=0.4;
PSO.dw=PSO.w1-PSO.w2;% init and final inertia weight
PSO.cg=2; PSO.cp=2;% social, congnitive factor
PSO.X=rand(PSO.D,PSO.N);
PSO.BND=[PSO.DLB PSO.DUB];% PSO bound

PSO.X(1,:)=PSO.BND(1,1)+PSO.X(1,:)*diff(PSO.BND(1,:)); % not sure what this does yet
history = zeros(1,PSO.G); % buffer for history

%% PSO main loop
carInit = car;
for g=1:PSO.G
    tic
    disp(PSO.X); % particles
    fit = zeros(1,PSO.N); % buffer for fit
    for n=1:PSO.N % number of particles
        car = carInit;
        d = PSO.X(n);  % pick one particle
        [virtual] = virtualNew(targetPosition,d);
        virtualVel = calculateDistance(virtual,target)/time.T;
        for t=0:time.dt:time.T
            vel = evalfis(calculateDistance(car,virtual),fis_vel);
            omega = evalfis(calculateAngle(car,virtual),fis_ang);
            [virtual] = CarNow(virtual,time,virtualVel,0);
            [car] = CarNow(car,time,vel,omega);
        end
        err = [car.x-target.x,car.y-target.y,wrapToPi(car.q-target.q)];
        fit(n) = sqrt(sum(err.^2));
        
        if fit(n) < PSO.gbest
            PSO.gbest = fit(n);
            PSO.Gbest = PSO.X(:,n);
        end
        if fit(n) < PSO.pbest(n)
            PSO.pbest(n) = fit(n);
            PSO.Pbest(:,n) = PSO.X(:,n);
        end
    end
    history(g) = PSO.gbest;
    fprintf('Generation %d Gbest %5.3f gbest %5.3f\n',g,PSO.Gbest.',PSO.gbest);% currently best results
    % PSO Update
    toc
    w=PSO.w2+(1-g/PSO.G)*PSO.dw;% decreasing inertia weight (w_initial-w_final) = PSO.dw
    PSO.V=w*rand(PSO.D,PSO.N).*PSO.V+...
    PSO.cp*rand(PSO.D,PSO.N).*(PSO.Pbest-PSO.X)+...
    PSO.cg*rand(PSO.D,PSO.N).*(repmat(PSO.Gbest,[1,PSO.N])-PSO.X);
    PSO.X=PSO.X+PSO.V;
    % repair particles
    z=find(PSO.X(1,:)<PSO.BND(1,1));
    PSO.X(1,z)=PSO.BND(1,1)+rand(1,length(z))*diff(PSO.BND(1,:));
    z=find(PSO.X(1,:)>PSO.BND(1,2));
    PSO.X(1,z)=PSO.BND(1,1)+rand(1,length(z))*diff(PSO.BND(1,:));
end

%% Showing the optimized distance
figure(fig);        % set to current axis
d = PSO.Gbest(1);   % get best distance
car = carInit;      % reset car position
% Generate Virtual target
virtualPosition = [targetPosition(1)-d*cos(targetHeading),targetPosition(2)-d*sin(targetHeading),targetHeading];
[virtual] = CarNew(virtualPosition,'k');
[virtual] = CarNow(virtual,time,0,0);

% plot all objects on fig
[virtual] = CarShow(fig,virtual,0);
[car] = CarShow(fig,car,0);
[target] = CarShow(fig,target,0);

%% Calculate the velocity at which the virtual target moves to real target
virtualVel = calculateDistance(virtual,target)/time.T;

for t=0:time.dt:time.T
    vel = evalfis(calculateDistance(car,virtual),fis_vel);
    omega = evalfis(calculateAngle(car,virtual),fis_ang);
    
    [virtual] = CarNow(virtual,time,virtualVel,0);
    [virtual] = CarShow(fig,virtual,t);
    
    [car]=CarNow(car,time,vel,omega);
    [car]=CarShow(fig,car,t);
end

%% Error vs Generation graph
figure();
plot(1:PSO.G,history);
xlabel('Generation');
ylabel('Error');
title('POS');
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
    car.q=wrapToPi(car.q);
    car.trace(:,end+1)=[car.x; car.y; car.q];
end

function [car]=CarShow(fig,car,t)
    ax=axis;
    Rz=[  cos(car.q) -sin(car.q); 
          sin(car.q)  cos(car.q)];
    shape=Rz*car.shape+repmat([car.x;car.y],1,6);
    set(car.hdL.shape,'xdata',shape(1,:),'ydata',shape(2,:)); 
    set(car.hdL.trace,'xdata',car.trace(1,:),'ydata',car.trace(2,:));
    axis(ax); title(sprintf('Time %d',t)); pause(0.001);
end

function [virtual]=virtualNew(target,d)
    % This function omits the plotting
    % to be used in POS loop
    virtual.x=target(1)-d*cos(target(3)); virtual.y=target(2)-d*sin(target(3)); virtual.q=target(3);
    virtual.trace=[virtual.x; virtual.y; virtual.q];
end

function distance = calculateDistance(src,dest)
    distance = sqrt((src.x-dest.x)^2+(src.y-dest.y)^2);
end

function omega = calculateAngle(src,dest)
    omega=atan2((dest.y-src.y),(dest.x-src.x))-src.q;
    omega=wrapToPi(omega);
end