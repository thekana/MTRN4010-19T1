
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

% load data file.
load(DataFileName); 
% The variable that contains the data structure is named "dataL"
% (because it was created under that name and saved)  

N = dataL.N;                         % number of scans in this squence.

disp('Brute force plotting');
disp('There is no GUI. Use "Control-C" to break the loop.');

for i=1:10:N                        % in this example, I skip some of them..
    scan_i = dataL.Scans(:,i);
    ProcessScan(scan_i);
    
    s=sprintf('(Brute force plotting...)\nShowing scan #[%d]/[%d]\r',i,N);
    title(s);
    
    %pause(0.01) ;                   % wait for ~10ms
end
disp('Done. Bye.');

return;
end


%.............................
function ProcessScan(scan)

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

% Plot. "BRUTE FORCE" plot (see note 1).
figure(1) ; clf(); hold on;
plot(X,Y,'b.');                     % all points
ii = find(intensities~=0);          % find those "pixels" that had intense reflection (>0) (aka: Highly Reflective pixels, HR)
plot(X(ii),Y(ii),'+r');             % plot highly reflective ones
axis([-10,10,-5,20]);                % focuses plot on this region ( of interest in L220)

dataForOOI = [X,Y,single(intensities)];
% To be done (by you)
OOIs = ExtractOOIs(dataForOOI);
PlotOOIs(OOIs);

return;
end

% find algo
function r = ExtractOOIs(dataForOOI)
    

    % your part...
    threshold = 0.15;
    distance = [0; sqrt(diff(dataForOOI(:,1)).^2+diff(dataForOOI(:,2)).^2)]; %vector of distance
    cluster_number = 1;
    cluster_vector = zeros(length(distance),1);
    for i=1:length(distance)
        cluster_vector(i)=cluster_number;
        if (distance(i) > threshold)
            cluster_number = cluster_number+1;
        end
    end
    %dataForOOI = [dataForOOI,cluster_vector];
    r.N = max(cluster_vector);
    r.Centers = zeros(2,r.N);
    r.Diameter = zeros(1,r.N);
    r.Color = zeros(1,r.N);
    for i = 1:r.N
        temp = dataForOOI(cluster_vector==i,:);
        [xc,yc,R] = circfit(temp(:,1),temp(:,2));
        r.Diameter(i) = R*2;
        r.Centers(1,i) = xc;
        r.Centers(2,i) = yc;
        r.Color(i) = max(temp(:,3))~=0;
    end
    
    Filter = r.Diameter >= 5 & r.Diameter <=20;
    
    r.Centers(:,~Filter)=[];
    r.Diameter(~Filter)=[];
    r.Color(~Filter)=[];
    r.N = length(r.Color);
    
return;
end
    
function PlotOOIs(OOIs)
    if OOIs.N<1, return ; end;
    % your part....
    plot(OOIs.Centers(1,:),OOIs.Centers(2,:),'c*');
    plot(OOIs.Centers(1,OOIs.Color>0),OOIs.Centers(2,OOIs.Color>0),'g*');
return;
end

function dist = getDistance(point1,point2)
    dist = sqrt((point2(1)-point1(1))^2+(point2(2)-point1(2))^2);
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
% --------------------------------
% note 1: for a more efficient way of dynamically plotting, see example
% "ExampleUseLaserData.m", where we used graphical objects handles.
% --------------------------------
% Questions?  Ask the lecturer, j.guivant@unsw.edu.au
% --------------------------------

