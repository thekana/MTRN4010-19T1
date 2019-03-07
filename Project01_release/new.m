
% Example program, for off-line processing of saved laser scans. 
% Example source code, useful for Task01. 
% AAS - 2018.S1.
% Jose Guivant.

% Note: Plotting of results is implemented, here, via a brute force approach.
% See other provided examples, for better implementation (dynamic updates using "set()").
% I expect you to avoid "brute force" implementations, in your projects. 


function MyProgram(DataFileName)

clc(); close all;
% In case the caller does not specify the input argument, we propose a default one.
if ~exist('DataFileName','var'), DataFileName ='Laser__2.mat'; end;

% Create a global variable struct to be used in multiple functions
global CCC;
CCC = []; CCC.flagPause = 0;

% load data file.
load(DataFileName); 
% The variable that contains the data structure is named "dataL"
% (because it was created under that name and saved)  

N = dataL.N;                        %number of scans in this squence.
figure('visible','on');
clf();
myHandle.handle1 = plot(0,0,'b.');  %handle for all points
hold on;
myHandle.handle2 = plot(0,0,'r+');  %handle for reflective points
axis([-10,10,0,20]);                %focuses plot on this region ( of interest in L220)
xlabel('x (meters)');
ylabel('y (meters)');
myHandle.handle3 = title('');
myHandle.handle4 = plot(0,0);       %handle for reflective OOIs
myHandle.handle5 = plot(0,0);       %handle for non-reflective OOIs
zoom on; grid on;
uicontrol('Style','pushbutton','String','Pause/Cont.','Position',[10,1,80,20],'Callback',{@PushButtonCallBack,1});

for i=1:N
    tic
    while (CCC.flagPause), pause(0.15); end
    scan_i = dataL.Scans(:,i);
    ProcessScan(scan_i,myHandle);
    
    s=sprintf('Showing scan #[%d]/[%d]\r',i,N);
    set(myHandle.handle3,'string',s);
    
    pause(0.01) ;                   % 10hz refresh rate
    toc
end
disp('Done. Bye.');

return;
end


%.............................
function ProcessScan(scan,myHandle)

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


set(myHandle.handle1,'xdata',X,'ydata',Y,'color','b','marker','.')         % update all points
ii = find(intensities~=0);          % find those "pixels" that had intense reflection (>0) (aka: Highly Reflective pixels, HR)
set(myHandle.handle2,'xdata',X(ii),'ydata',Y(ii),'color','r','marker','+');     % plot highly reflective ones

data = [X,Y,single(intensities)];   % concat vectors

OOIs = ExtractOOIs(data);
PlotOOIs(OOIs,myHandle);

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
    
%     for i = 1:r.N
%         cluster_i = data(cluster_vector==i,:);
%         r.Centers(1,i) = mean(cluster_i(:,1));
%         r.Centers(2,i) = mean(cluster_i(:,2));
%         r.Color(i) = max(cluster_i(:,3))~=0;
%         [m,~] = size(cluster_i);
%         if m == 1
%             r.Diameter(i) = 0;
%         else     
%             r.Diameter(i) = max(pdist(cluster_i(:,1:2)));
%         end
%     end
    % Filling r struct using circfit function from Izhak bucher 25/oct /1991
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
    set(myHandle.handle4,'xdata',OOIs.Centers(1,OOIs.Color>0),'ydata',OOIs.Centers(2,OOIs.Color>0),'color','g','marker','+','markersize',10);
    set(myHandle.handle5,'xdata',OOIs.Centers(1,OOIs.Color==0),'ydata',OOIs.Centers(2,OOIs.Color==0),'color','k','marker','+','markersize',10);
return;
end

function   [xc,yc,R] = circfit(x,y)
%
%   [xc yx R] = circfit(x,y)
%
%   fits a circle  in x,y plane in a more accurate
%   (less prone to ill condition )
%  procedure than circfit2 but using more memory
%  x,y are column vector where (x(i),y(i)) is a measured point
%
%  result is center point (yc,xc) and radius R
%  an optional output is the vector of coeficient a
% describing the circle's equation
%
%   x^2+y^2+a(1)*x+a(2)*y+a(3)=0
%
%  By:  Izhak bucher 25/oct /1991, 
   x=x(:); y=y(:);
   a=[x y ones(size(x))]\[-(x.^2+y.^2)];
   xc = -.5*a(1);
   yc = -.5*a(2);
   R  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
end


function PushButtonCallBack(~,~,x)
    global CCC;  
    if (x==1)
       CCC.flagPause = ~CCC.flagPause; %Switch ON->OFF->ON -> and so on.
       disp(x);
       disp(CCC.flagPause);
    end
end