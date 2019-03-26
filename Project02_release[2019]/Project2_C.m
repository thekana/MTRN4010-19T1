
clc(); close all;


% Create a global variable struct to be used in multiple functions
global CCC;
CCC = []; CCC.flagPause = 0;
global landmark;
landmark.coor = [];
landmark.id = [];
% load data file.
load('DataForProject02/IMU_dataC.mat');
load('DataForProject02/Speed_dataC.mat');
load('DataForProject02/Laser__2C.mat');

% The variable that contains the data structure is named "dataL"
% (because it was created under that name and saved)

% Obtain an estimate of the position and orientation of the platform at
% at scan's timestamp
[x,y,theta] = GetData(IMU,Vel,dataL);


N = dataL.N;                        %number of scans in this squence.
figure('visible','on');
clf();
%myHandle.handle1 = plot(0,0,'b.');  %handle for all points
hold on;
%myHandle.handle2 = plot(0,0,'r+');  %handle for reflective points
axis([-10,10,0,20]);                %focuses plot on this region ( of interest in L220)
xlabel('x (meters)');
ylabel('y (meters)');
myHandle.handle3 = title('');
myHandle.handle4 = plot(0,0);       %handle for reflective OOIs
myHandle.handle5 = plot(0,0);       %handle for non-reflective OOIs
landmark.handle = plot(0,0,'linestyle','none');
zoom on; grid on;
uicontrol('Style','pushbutton','String','Pause/Cont.','Position',[10,1,80,20],'Callback',{@PushButtonCallBack,1});

for i=1:N
    
    while (CCC.flagPause), pause(0.15); end
    scan_i = dataL.Scans(:,i);
    tic
    ProcessScan(scan_i,myHandle,x(i),y(i),theta(i));
    toc
    s=sprintf('Showing scan #[%d]/[%d]\r',i,N);
    set(myHandle.handle3,'string',s);
    
    pause(0.01) ;                   % 10hz refresh rate
    
end
disp('Done. Bye.');


function [xKL,yKL,thetaKL] = GetData(IMU,Vel,dataL)
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
    % thetaK = thetaK * 180/pi;
    % thetaKL = thetaKL * 180/pi;
end

%.............................
function ProcessScan(scan,myHandle,x,y,theta)
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
    OOIs = ToGlobalCoordinateFrame(OOIs,x,y,theta);
    PlotOOIs(OOIs,myHandle);
    IdentifyOOIs(OOIs);

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
    Filter = (r.Diameter >= 0.05 & r.Diameter <=0.20);
    %Clearing non OOI objects from r then resize r
    r.Centers(:,~Filter)=[];
    r.Diameter(~Filter)=[];
    r.Color(~Filter)=[];
    r.N = length(r.Color);
    
return;
end

function PlotOOIs(OOIs,myHandle)
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
    
    d = 0.46;
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

function IdentifyOOIs(r)
    global landmark;
    OOIarray = r.Centers(:,r.Color>0);
    [~,n] = size(OOIarray);
    if isempty(landmark.coor)
        % add all ooi and give unique id
        landmark.coor = OOIarray;
        landmark.id = 1:length(landmark.coor);
    elseif ~isempty(OOIarray)
        % matching process
        for i = 1:n 
            match = 0;
            for j = landmark.id
                distance = norm(OOIarray(:,i)-landmark.coor(:,j));
                if (distance <= 0.4)
                    match = 1;
                end
            end
            if ~match
                landmark.coor = [landmark.coor,OOIarray(:,i)];
                landmark.id = [landmark.id,j+1];
            end
        end    
    end     
    set(landmark.handle,'xdata',landmark.coor(1,:),'ydata',landmark.coor(2,:),'color','k','marker','+','markersize',10);
end

    