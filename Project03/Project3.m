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
axis([-5,5,0,10]);                %focuses plot on this region ( of interest in L220)
xlabel('x (meters)');
ylabel('y (meters)');
global myHandle;
myHandle.handle3 = title('');
myHandle.handle4 = plot(0,0);       %handle for reflective OOIs
myHandle.handle5 = plot(0,0);       %handle for non-reflective OOIs
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


%%Populating Xe using EKF
DataHandling(IMU,Vel,dataL);


%%For loop to plot I guess
% for i=1:N
%     
%     while (CCC.flagPause), pause(0.15); end
%     s=sprintf('Showing scan #[%d]/[%d]\r',i,N);
%     set(myHandle.handle3,'string',s);
%     plotRobot(robot,x(i),y(i),theta(i));
%     
%     pause(0.01) ;                   % 10hz refresh rate
%     
% end
disp('Done. Bye.');


function DataHandling(IMU,Vel,dataL)
    %%% Noises
    stdDevGyro = 1.4;
    stdDevSpeed = 0.4;
    stdRangeMeasure = 0.16;
    stdBearingMeasure = 1.1;
    
    global landmark;
    global myHandle;
    global CCC;
    %%Matrices
    P = zeros(4,4);
    P(4,4) = (4*pi/180)^2;
    Pu = diag([stdDevSpeed^2,stdDevGyro^2]);
    Q1 = diag([ (0.01)^2 ,(0.01)^2 , (1*pi/180)^2,0]);
    %%
    % remove bias
    time = double(IMU.times-IMU.times(1))/10000;
    Laser_time = double(dataL.times-dataL.times(1))/10000;

    yaw = -IMU.DATAf(6,:); % negate w(K)

    horizon0 = length(find(time<20)); %stationary for 20 secs
    Bias = mean(yaw(1:horizon0)); %Bias

    yawC = yaw - Bias; % w(K) corrected
    L = length(yawC);

    % buffer for variables temp 40K length
    Xtemp = zeros(4,L);
    Xtemp(:,1) =  [0;0;pi/2;0];
    % buffer for what we actually care about
    Xe = zeros(4,length(Laser_time));
    Xe(:,1) = [0;0;pi/2;0];
    
    j = 2;
    
    % must obtain the robot position and heading at the time scan[i]
    for i = 2:L
        dt = time(i)-time(i-1); % find dt
        Xtemp(3,i) = Xtemp(3,i-1) + dt*yawC(i-1); 
        Xtemp(1:2,i) = Xtemp(1:2,i-1) + dt*[Vel.speeds(i-1)*cos(Xtemp(3,i));Vel.speeds(i-1)*sin(Xtemp(3,i))];
        Xtemp(4,i) = Xtemp(4,i-1) + dt*0;
        % handling laser scan and position frequency difference
        if (j <= length(Laser_time) && Laser_time(j) - time(i-1)< dt)
            dtL = Laser_time(j) - time(i-1);
            %%Matrices
            J = [ [1,0,-dtL*Vel.speeds(i-1)*sin(Xe(3,j-1)),0]  ; [0,1,dtL*Vel.speeds(i-1)*cos(Xe(3,j-1)),0];[ 0,0,1,-dtL]; [0,0,0,1]]; %4x4
            Ju = [dtL*cos(Xe(3,j-1)),0;dtL*sin(Xe(3,j-1)),0;0,dtL;0,0]; %4x2
            Q = Ju*Pu*Ju'+Q1;
            P = J*P*J'+Q ;
            %%
            Xe(:,j) = Xtemp(:,i-1) + dtL*[Vel.speeds(i-1)*cos(Xtemp(3,i));Vel.speeds(i-1)*sin(Xtemp(3,i)); yawC(i-1);0];
            %ProcessScan(dataL.Scans(:,j-1),Xe(1,j-1),Xe(2,j-1),Xe(3,j-1));
            Local_OOIs = ProcessScan(dataL.Scans(:,j-1));
            Global_OOIs = ToGlobalCoordinateFrame(Local_OOIs,Xe(1,j-1),Xe(2,j-1),Xe(3,j-1)); % convert OOIs to global coor
            PlotOOIs(Global_OOIs); %Will plot only bright points
            DataAssociation(Global_OOIs,Local_OOIs); % identify landmark && data association
            if(landmark.detected >0)
                %EKF parts
                for u = 1:landmark.detected
                    ID = landmark.DAcoor(3,u);
                    d = 0.46;
                    %Calc expected
                    eDX = (landmark.coor(1,ID)-Xe(1,j)+d*cos(Xe(3,j)));
                    eDY = (landmark.coor(2,ID)-Xe(2,j)+d*sin(Xe(3,j)));
                    eDD = sqrt( eDX*eDX + eDY*eDY );
                    H = [  -eDX/eDD , -eDY/eDD , 0,0;
                        eDY/eDD^2, -eDX/eDD^2, -1,0]; 
                    ExpectedRange = eDD;
                    ExpectedAngle = atan2(eDY,eDX) - Xe(3,j) + pi/2;
                    %Calc measurement
%                     eMX = (landmark.DAcoor(1,u)-Xe(1,j));
%                     eMY = (landmark.DAcoor(2,u)-Xe(2,j));
%                     eMD = sqrt( eMX*eMX + eMY*eMY );
%                     MeasuredRange = eMD;
%                     MeasuredAngle = atan2(eMY,eMX) - Xe(3,j) + pi/2;

                    eMX = (landmark.localOOIs(1,u));
                    eMY = (landmark.localOOIs(2,u));
                    eMD = sqrt( eMX*eMX + eMY*eMY );
                    MeasuredRange = eMD;
                    MeasuredAngle = atan2(eMY,eMX) - Xe(3,j) + pi/2;
                    z = [MeasuredRange-ExpectedRange;
                        wrapToPi(MeasuredAngle-ExpectedAngle)];
                    R = diag([stdRangeMeasure^2*4 stdBearingMeasure^2*4]);
                    %%EKFSteps
                    S = R + H*P*H';
                    iS=inv(S);
                    K = P*H'*iS;
                    Xe(:,j) = Xe(:,j)+K*z;
                    P = P-K*H*P;
                end    
            end
            while (CCC.flagPause), pause(0.15); end
            s=sprintf('Showing scan #[%d]/[%d]\r',j,length(Laser_time));
            set(myHandle.handle3,'string',s);
            plotRobot(Xe(1,j-1),Xe(2,j-1),Xe(3,j-1));
            pause(0.01) ;                   % 10hz refresh rate
            j = j + 1;
        end
    end
        set(myHandle.real,'xdata',Xtemp(1,:),'ydata',Xtemp(2,:),'color','r')
    % convert from radian to degree
    % thetaK = thetaK * 180/pi;
    % thetaKL = thetaKL * 180/pi;
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
            if (distance <= 0.4)
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
    set(robot.heading,'xdata',coor(1,:),'ydata',coor(2,:),'markersize',2,'color','y');
    set(robot.trace,'xdata',robot.traceData(1,:),'ydata',robot.traceData(2,:),'color','b');
end
    