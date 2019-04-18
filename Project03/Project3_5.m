clc(); close all;
%% Create a global variable struct to be used in multiple functions
global CCC;
CCC = []; CCC.flagPause = 0;
global landmark;
landmark.coor = [];
landmark.id = [];
landmark.detected = 0;
landmark.DAcoor = [];
landmark.localOOIs = [];
%% Setup loading data file.
load('DataForProject02/IMU_dataC.mat');
load('DataForProject02/Speed_dataC.mat');
load('DataForProject02/Laser__2C.mat');

% The variable that contains the data structure is named "dataL"
% (because it was created under that name and saved)


N = dataL.N;                        %number of scans in this squence.
figure('visible','on');
clf();
hold on;
axis([-5,3,-1,7]);                %focuses plot on this region ( of interest in L220)
xlabel('x (meters)');
ylabel('y (meters)');
global myHandle;
myHandle.handle3 = title('');
myHandle.handle4 = plot(0,0);       %handle for reflective OOIs
myHandle.handle5 = plot(0,0);       %handle for non-reflective OOIs
myHandle.handle6 = plot(0,0);       %for DR
myHandle.real = plot(0,0);
%% Creating landmark and currently detected OOIs graphic handles
landmark.handle = plot(0,0,'linestyle','none');
landmark.text = text(0,0,'');
landmark.DA = text(zeros(1,5),zeros(1,5),'');
%% Creating robot body
global robot;

robot.trace = plot(0,0);
robot.traceData = [];
robot.body = plot(0,0);
robot.heading = plot(0,0);
%%
zoom on; grid on;
uicontrol('Style','pushbutton','String','Pause/Cont.','Position',[10,1,80,20],'Callback',{@PushButtonCallBack,1});

%%% Noises
stdDevGyro = 1.4;
stdDevSpeed = 0.4;
stdRangeMeasure = 0.16;
stdBearingMeasure = 1.1;
max_acceleration = 1.5;

b = 100*pi/180; %between -2,2
%%Matrices
P = zeros(5,5);
P(4,4) = b^2;
Pu = diag(stdDevGyro^2);

time = double(IMU.times-IMU.times(1))/10000;
Laser_time = double(dataL.times-dataL.times(1))/10000;

yaw = -IMU.DATAf(6,:); % negate w(K)

horizon0 = length(find(time<20)); %stationary for 20 secs
Bias = mean(yaw(1:horizon0)); %Bias

yawC = yaw;
L = length(yawC);

% buffer for variables temp 40K length
Xehistory = zeros(5,L);
Xdrhistory = zeros(3,L);
Xe =  [0;0;pi/2;0;0];
Xdr = [0;0;pi/2];
% buffer for what we actually care about
Ltime = length(Laser_time);

current_scan = 1;

for i = 2:length(time)-1
    dt = time(i)-time(i-1); % find dt
    imuGyro = yawC(i)-Bias;
    speed = Vel.speeds(i);
    Xdr = processModelDR(imuGyro,speed,dt,Xdr);
    Xdrhistory(:,i) = Xdr;
end
% must obtain the robot position and heading at the time scan[i]
for i = 2:length(time)-1
    dt = time(i)-time(i-1); % find dt
    imuGyro = yawC(i);
    detectedOOIs = 0;
    Q1 = diag([ (0.01)^2 ,(0.01)^2 , (1*pi/180)^2,0,(dt*max_acceleration)^2]); %5x
    %% Process scan when there is laser data
    if (current_scan <= length(Laser_time) && Laser_time(current_scan) - time(i)< dt)
        % Get X,Y to process scans and transform OOIs
        dtL = Laser_time(current_scan) - time(i-1);
        
        Xtemp = processModel(imuGyro,dtL,Xe);
        
        Local_OOIs = ProcessScan(dataL.Scans(:,current_scan));
        Global_OOIs = ToGlobalCoordinateFrame(Local_OOIs,Xtemp(1),Xtemp(2),Xtemp(3)); % convert OOIs to global coor
        PlotOOIs(Global_OOIs); %Will plot only bright points
        DataAssociation(Global_OOIs,Local_OOIs); % identify landmark && data association
        detectedOOIs = landmark.detected;
        current_scan = current_scan + 1;
    end
    
    J = [ [1,0,-dt*Xe(5)*sin(Xe(3)),0,dt*cos(Xe(3))] ; [0,1,dt*Xe(5)*cos(Xe(3)),0,dt*sin(Xe(3))];[ 0,0,1,-dt,0];[0,0,0,1,0];[0,0,0,0,1]]; %5x5 jacobian
    Ju = [0;0;dt;0;0]; %5x1 linear transformation of the input which is now only gyro
    Q = Ju*Pu*Ju'+Q1;
    P = J*P*J'+Q ;
    Xe = processModel(imuGyro,dt,Xe);
    
    if(detectedOOIs >0)
        %EKF parts
        for u = 1:landmark.detected
            ID = landmark.DAcoor(3,u);
            d = 0.46;
            %Calc expected
            eDX = (landmark.coor(1,ID)-Xe(1)-d*cos(Xe(3)));
            eDY = (landmark.coor(2,ID)-Xe(2)-d*sin(Xe(3)));
            eDD = sqrt( eDX*eDX + eDY*eDY );
            H = [  -eDX/eDD , -eDY/eDD , 0,0,0; %2x5
                eDY/eDD^2, -eDX/eDD^2, -1,0,0]; 
            ExpectedRange = eDD;
            ExpectedAngle = atan2(eDY,eDX) - Xe(3) + pi/2;
            
            eMX = (landmark.localOOIs(1,u));
            eMY = (landmark.localOOIs(2,u));
            eMD = sqrt( eMX*eMX + eMY*eMY );
            MeasuredRange = eMD;
            MeasuredAngle = atan2(eMY,eMX); 
            z = [MeasuredRange-ExpectedRange;
                wrapToPi(MeasuredAngle-ExpectedAngle)];
            R = diag([stdRangeMeasure^2*4 stdBearingMeasure^2*4]);
            %%EKFSteps
            S = R + H*P*H';
            iS=inv(S);
            K = P*H'*iS;
            Xe = Xe+K*z;
            P = P-K*H*P;

        end
    end
    Xehistory(:,i) = Xe;
    while (CCC.flagPause), pause(0.15); end
    s=sprintf('Showing scan #[%d]/[%d]\tEstimated Bias[%f] Velocity[%f]\r',i+1,length(time),Xe(4),Xe(5));
    set(myHandle.handle3,'string',s);
    set(myHandle.handle6,'xdata',Xdrhistory(1,1:i),'ydata',Xdrhistory(2,1:i),'LineStyle','none','marker','.');
    plotRobot(Xe(1),Xe(2),Xe(3));
    pause(0.00001) ;                   % 10hz refresh rate
end
    figure()
    plot(time(1:length(time)-1),Xehistory(4,1:length(time)-1));
    axis([0,250,0,0.02]);
    figure()
    plot(time(1:length(time)-1),Xehistory(5,1:length(time)-1));
    hold on
    plot(time,Vel.speeds)
function Xnext = processModel(omega,dt,Xprev)
    
    Xnext = zeros(5,1); 
    Xnext(1) = Xprev(1) + Xprev(5)*cos(Xprev(3))*dt;
    Xnext(2) = Xprev(2) + Xprev(5)*sin(Xprev(3))*dt;
    Xnext(3) = Xprev(3) + dt*(omega-Xprev(4));
    Xnext(4) = Xprev(4) + 0;
    Xnext(5) = Xprev(5) + 0;
end

function Xnext = processModelDR(omega,speed,dt,Xprev)
    
    Xnext = zeros(3,1); 
    Xnext(1) = Xprev(1) + speed*cos(Xprev(3))*dt;
    Xnext(2) = Xprev(2) + speed*sin(Xprev(3))*dt;
    Xnext(3) = Xprev(3) + omega*dt;  
end

%.............................
function OOIs=ProcessScan(scan)
    % Extract range and intensity information, from raw measurements.
    % Each "pixel" is represented by its range and intensity of reflection.
    % It is a 16 bits number whose bits 0-12 define the distance (i.e. the range)
    % in cm (a number 0<=r<2^13), and bits 13-15 indicate the intensity 
    %( a number 0<=i<8 ).
    % We extract range and intensity, here.
    %useful masks, for dealing with bits.
    mask1FFF = uint16(2^13-1);
    maskE000 = bitshift(uint16(7),13)  ;

    intensities = bitand(scan,maskE000);

    ranges    = single(bitand(scan,mask1FFF))*0.01; 
    % Ranges expressed in meters, and unsing floating point format (e.g. single).

    % 2D points, expressed in Cartesian. From the sensor's perpective.
    angles = [0:360]'*0.5* pi/180 ;         % associated angle, for each individual range in a scan
    X = cos(angles).*ranges;
    Y = sin(angles).*ranges;    


    %set(myHandle.handle1,'xdata',X,'ydata',Y,'color','b','marker','.')         % update all points
    % ii = find(intensities~=0);          % find those "pixels" that had intense reflection (>0) (aka: Highly Reflective pixels, HR)
    %set(myHandle.handle2,'xdata',X(ii),'ydata',Y(ii),'color','r','marker','+');     % plot highly reflective ones

    data = [X,Y,single(intensities)];   % concat vectors

    OOIs = ExtractOOIs(data);

return;
end

function r = ExtractOOIs(data)
   
    threshold = 0.15;
    distance = [0; sqrt(diff(data(:,1)).^2+diff(data(:,2)).^2)]; %vector of distance between two consecutive points
    cluster_number = 1; % starting cluster number
    cluster_vector = zeros(length(distance),1);
    
    for i=1:length(distance)
        % assigning cluster number to points
        if distance(i) >= threshold
            cluster_number = cluster_number+1;
        end
        cluster_vector(i)=cluster_number;
    end
    r.N = max(cluster_vector);
    r.Centers = zeros(2,r.N);
    r.Diameter = zeros(1,r.N);
    r.Color = zeros(1,r.N);
    % Filling r struct using circfit function from Izhak bucher 25/oct /1991
    % https://au.mathworks.com/matlabcentral/fileexchange/5557-circle-fit
    for i = 1:r.N
        cluster_i = data(cluster_vector==i,:);
        % checking matrix size discard every cluster that has fewer than 3
        % points
        size_check = size(cluster_i);
        if size_check(1) < 3
           %disp('skip');
           continue;
        end
        [xc,yc,R] = circfit(cluster_i(:,1),cluster_i(:,2));
        r.Centers(1,i) = xc;
        r.Centers(2,i) = yc;
        r.Diameter(i) = 2*R;
        r.Color(i) = max(cluster_i(:,3))~=0;
    end
    Filter = (r.Diameter >= 0.03 & r.Diameter <=0.25);
    %Clearing non OOI objects from r then resize r
    r.Centers(:,~Filter)=[];
    r.Diameter(~Filter)=[];
    r.Color(~Filter)=[];
    r.N = length(r.Color);
    
return;
end

function PlotOOIs(OOIs)
    global myHandle;
    if OOIs.N<1, return ; end;
    myHandle.handle4.LineStyle = 'none';
    myHandle.handle4.LineWidth = 2;
    myHandle.handle5.LineStyle = 'none';
    myHandle.handle5.LineWidth = 1.5;
    set(myHandle.handle4,'xdata',OOIs.Centers(1,OOIs.Color>0),'ydata',OOIs.Centers(2,OOIs.Color>0),'color','g','marker','*','markersize',10);
    %set(myHandle.handle5,'xdata',OOIs.Centers(1,OOIs.Color==0),'ydata',OOIs.Centers(2,OOIs.Color==0),'color','k','marker','+','markersize',10);
return;
end

function PushButtonCallBack(~,~,x)
    global CCC;  
    if (x==1)
       CCC.flagPause = ~CCC.flagPause; %Switch ON->OFF->ON -> and so on.
       disp(x);
       disp(CCC.flagPause);
    end
end

function OOIs = ToGlobalCoordinateFrame(OOIs,x,y,theta)
    
    d = 0.46; % from project specs
    % negate x-coordinate
    OOIs.Centers(1,:) = -OOIs.Centers(1,:);
    % adding d
    OOIs.Centers(2,:) = d + OOIs.Centers(2,:);
    a = theta - pi/2;
    R = [ cos(a), -sin(a);
          sin(a), cos(a)];
    %OOI
    OOIs.Centers = R * OOIs.Centers + [x;y];
    
end

function DataAssociation(r,rLocal)
    global landmark;
    OOIglobal = r.Centers(:,r.Color>0);
    OOIlocal = rLocal.Centers(:,rLocal.Color>0);
    %negate Xaxis in local OOI
    OOIlocal(1,:) = -OOIlocal(1,:);
    [~,n] = size(OOIglobal);
    DA = []; %data association
    Measured = [];
    landmark.detected = 0;
    if isempty(landmark.coor)
        % add all ooi and give unique id
        landmark.coor = OOIglobal;
        landmark.id = 1:length(landmark.coor);
        landmark.text = text(landmark.coor(1,:)+0.1,landmark.coor(2,:)-0.1,string(landmark.id));
    end
    for i = 1:n 
        for j = landmark.id
            distance = norm(OOIglobal(:,i)-landmark.coor(:,j));
            if (distance <= 0.8)
                temp = [OOIglobal(:,i);j];
                DA = [DA,temp];
                temp = [OOIlocal(:,i);j];
                Measured = [Measured,temp];
            end
        end
        landmark.DAcoor = DA;
        landmark.localOOIs = Measured;
        if ~isempty(DA)
            landmark.detected = length(DA(3,:));
        end
    end    
    
    set(landmark.handle,'xdata',landmark.coor(1,:),'ydata',landmark.coor(2,:),'color','k','marker','+','markersize',10);
    
    for i = 1:5
        set(landmark.DA(i),'String','');
        if ~isempty(DA)
           if i<=length(DA(3,:))
               set(landmark.DA(i),'Position',DA(1:2,i),'String',string(DA(3,i)));
           end
        end    
    end
end

function plotRobot(x,y,theta)
    global robot;
    robot.traceData = [robot.traceData,[x;y]];
    set(robot.body,'xdata',x,'ydata',y,'markersize',7,'marker','diamond');
    R = [ cos(theta), -sin(theta);
          sin(theta), cos(theta)];
    coor = [0 0.2 0.4 0.8 1; 0 0 0 0 0];
    coor = R * coor + [x;y];
    set(robot.heading,'xdata',coor(1,:),'ydata',coor(2,:),'markersize',2,'color','k');
    set(robot.trace,'xdata',robot.traceData(1,:),'ydata',robot.traceData(2,:),'color','b');
end
    