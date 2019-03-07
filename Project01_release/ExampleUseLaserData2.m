
% AAS-mtrn4010. 2019.S1.  Useful code for Project01.
% (modified version, which includes some basic control buttons)

% Example program, for processing laser scans, stored in a Matlab data file.
% It shows the scans in the original/native way :POLAR.
% You are requested to modify this program to show the data in Cartesian.


% The data correspond to a laser scanner  (LIDAR) installed at the front of our UGV
% (robot). It is pointing ahead the robot, scanning horizontally, i.e. in 2D
% When you plot the data in Cartesian you will see how the room did look
% from the perspective of the moving platform.


% IMPORTANT ==> Read the program and comments, before trying to modify it, and before asking. 


function ExampleUseLaserData2(file)

global ABCD;
ABCD.flagPause=0;
% In case the caller does not specify the input argument, we propose a
% default one, assumed to be in the same folder where we run this example from.
if ~exist('file','var'), file ='Laser__2.mat'; end;

load(file); 
% now, after loading, the data is in a variable named "dataL";



    % --------------------------------------
    % Create graphical object for refreshing data during program execution.
    figure(1) ; clf(); 
    
    MyGUIHandles.handle1 = plot(0,0,'b.');      % to be used for showing the laser points
   
    
    axis([0,180,0,20]);                         % focuses the plot on this region (of interest, close to the robot)
    xlabel('angle (degrees)');
    ylabel('range (meters)');
    
    MyGUIHandles.handle2 = title('');           % create an empty title..
    zoom on ;  grid on;
    
    % If you do not understand these functions ( figure() ,plot(),
    % axis(),.....) ===> Read Matlab's Help.
    % (In MTRN2500 you used them)
    
    %---------------------------------
    
    disp('Showing laser scans, in POLAR representation');
    disp('Yo need to modify this program for showing the data in Cartesian.');
    disp('(then the images will make sense, for our brains)');
    
    fprintf('\nThere are [ %d ] laser scans in this dataset (file [%s])\n',dataL.N,file);
    
    
    uicontrol('Style','pushbutton','String','Pause/Cont.','Position',[10,1,80,20],'Callback',{@MyCallBackA,1});
    uicontrol('Style','pushbutton','String','END Now','Position',[90,1,80,20],'Callback',{@MyCallBackA,2});

    
    
    % Now, loop through the avaialable scans..
    
N = dataL.N; 
skip=3;
i=1;
while 1,             % in this example I skip some of the laser scans.
    
    if (ABCD.flagPause), pause(0.2) ; continue ; end;
    if i>N, break ;  end;
    
        
        
    % Native time expressed via uint32 numbers, where 1 unint means
    % [1/10,000]second (i.e. 0.1 millisecond)
    % (BTW: we do not use the time, in this task)
    t =  double(dataL.times(i)-dataL.times(1))/10000;
    % t: time expressed in seconds, relative to the time of the first scan.
    
    scan_i = dataL.Scans(:,i);
    MyProcessingOfScan(scan_i,t,MyGUIHandles,i);   % some function to use the data...

    pause(0.01) ;                   % wait for ~10ms (approx.)
    i=i+skip;
end;

fprintf('\nDONE!\n');

% Read Matlab Help for explanation of FOR loops, and function double( ) and pause()


return;
end
%-----------------------------------------
function MyProcessingOfScan(scan,t,mh,i)
    % I made this function, to receive the following parameters/variables:
    % 'scan' : scan measurements to be shown.
    % 't':  associated time.    
    % 'i' : scan number.
    % 'mh'  : struct contaning handles of necessary graphical objects.
    
    angles = [0:360]'*0.5 ;         % Associated angle for each range of scan
    % same as in "dataL.angles".
    
    % scan data is provided as a array of class uint16, which encodes range
    % and intensity (that is the way the sensor provides the data, in that
    % mode of operation)
     
    MaskLow13Bits = uint16(2^13-1); % mask for extracting the range bits.
    % the lower 13 bits are for indicating the range data (as provided by this sensor)
    
    rangesA = bitand(scan,MaskLow13Bits) ; 
    % rangesA now contains the range data of the scan, expressed in CM, having uint16 format.
    
    % now I convert ranges to meters, and also to floating point format
    ranges    = 0.01*double(rangesA); 

    
    % and then refresh the data in the plots...
    set(mh.handle1,'xdata',angles,'ydata',ranges);

    % and some text...
    s= sprintf('Laser scan # [%d] at time [%.3f] secs',i,t);
    set(mh.handle2,'string',s);
    
    
    % Use Matlab help for learning the functionality of  uint16(), bitand()
    % , set( ), sprintf()...
    
    return;
end
% ---------------------------------------
% Callback function. I defined it, and associated it to certain GUI button,
function MyCallBackA(~,~,x)   
    global ABCD;
        
    if (x==1)
       ABCD.flagPause = ~ABCD.flagPause; %Switch ON->OFF->ON -> and so on.
       return;
    end;
    if (x==2)
        
        disp('you pressed "END NOW"');
        uiwait(msgbox('Ooops, you still need to implement this one!','?','modal'));
        
        % students complete this.
        return;
    end;
    return;    
end



%-----------------------------------------
    
% The data in the file is stored according to certain structure. 
% The fields we need are: "Scans" and "times".
% The number of scans is indicated by the field "N".
% The scans themselves are stored in the field "Scans", that is a matrix of 361 x N 
% (because each scan contains 361 measurements, that cover 180 degrees at steps of 1/2 degree).
% the ranges are originally expressed in Centimeters, as 13 bits unsigned int integers.
% as explained in the code.

%-----------------------------------------% 

%  by Jose Guivant.  MTRN4010.S1.2019 

%-----------------------------------------% 

%  Questions :Moodle or email to the lecturer.

%-----------------------------------------
