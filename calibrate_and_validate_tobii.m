%% Calibration and validation of Tobii Eyetracker in Matlab 
% Needs psychtoolbox (aka PTB) and gstreamer installed
% Tested on:
%   Matlab R2017b
%   PsychtoolboxVersion '3.0.15 - Flavor: beta - Corresponds to SVN Revision 9685
%   gstreamer-1.0-mingw-x86_64-1.16.1

% NOTE 1: Psychtoolbox takes f******g ages to start after windows restart...
% patience !
%
% Note 2: Few changes were made in SDKs provided by Tobii-> do not replace the
% SDK's folder in this directiory

% Based on:
% http://developer.tobiipro.com/matlab/matlab-sdk-reference-guide.html
% https://github.com/tobiipro/prosdk-addons-matlab
% 
% Coded by Adam 2019 
 
clear;
clc;
%% Add Tobii SDK to path
addpath(genpath('TobiiPro.SDK.Matlab_1.7.1.4'));
addpath('ScreenBasedCalibrationValidation');
 
%% Find All EyeTrackers
try
    disp(' ')
    disp('Connecting to Tobii:')
    disp('-------------------')

    max_connection_attempts=0;
    i=0; 
    while i<=max_connection_attempts
        Tobii = EyeTrackingOperations();
        eyetrackers = Tobii.find_all_eyetrackers();
        for i=1:size(eyetrackers,2)
            disp('Device found:')
            disp(['Address:',eyetrackers(i).Address])
            disp(['Name:',eyetrackers(i).Name])
            disp(['Serial Number:',eyetrackers(i).SerialNumber])
        end
        if ~isempty(eyetrackers); break; end % stop looking for tobii when found
        i=i+1;
        pause(0.1)
    end
    eyetracker=eyetrackers(1);    
catch
    disp('Tobii eyetracker cannot be found. Is it ON and connected?')
    error('Cannot connect to Tobii Eyetracker. Is it ON and connected?')
end

%% Set up PsychoToolbox
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);
Screen('Preference','SuppressAllWarnings',1 ); % this disables all messages from PTB; uncooment when something is broken :)
screens = Screen('Screens');

% So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% this. For help see: help max
screenNumber = max(screens);

HideCursor(screenNumber); % hides mouse cursor- looks a bit better

% Define black and white (white will be 1 and black 0). This is because
% luminace values are (in general) defined between 0 and 1.
% For help see: help WhiteIndex and help BlackIndex
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);

% Open an on screen window and color it black.
% For help see: Screen OpenWindow?
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Get the size of the on screen window in pixels.
% For help see: Screen WindowSize?
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
screen_pixels = [screenXpixels  screenYpixels];

% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[xCenter, yCenter] = RectCenter(windowRect);

% Enable alpha blending for anti-aliasing
% For help see: Screen BlendFunction?
% Also see: Chapter 6 of the OpenGL programming guide
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

%% Show tracking area and check user's distance from the screen
% Show on screen the position of the user's eyes to garantee that the user is positioned correctly in the tracking area.

% Dot size in pixels
dotSizePix = 30;

% Start collecting data
% The subsequent calls return the current values in the stream buffer.
% If a flat structure is prefered just use an extra input 'flat'.
% i.e. gaze_data = eyetracker.get_gaze_data('flat');
eyetracker.get_gaze_data();

Screen('TextSize', window, 20);

while ~KbCheck
    
    DrawFormattedText(window, 'When correctly positioned press any key to continue.', 'center', screenYpixels * 0.1, white);
    
    distance = [];
    
    gaze_data = eyetracker.get_gaze_data();
    
    if ~isempty(gaze_data)
        last_gaze = gaze_data(end);
        
        validityColor = [1 0 0];
        
        % Check if user has both eyes inside a reasonable tacking area.
        if last_gaze.LeftEye.GazeOrigin.Validity.value && last_gaze.RightEye.GazeOrigin.Validity.value %EDITED ADAM last_gaze.LeftEye.GazeOrigin.Validity && last_gaze.RightEye.GazeOrigin.Validity
            left_validity = all(last_gaze.LeftEye.GazeOrigin.InTrackBoxCoordinateSystem(1:2) < 0.85) ...
                && all(last_gaze.LeftEye.GazeOrigin.InTrackBoxCoordinateSystem(1:2) > 0.15);
            right_validity = all(last_gaze.RightEye.GazeOrigin.InTrackBoxCoordinateSystem(1:2) < 0.85) ...
                && all(last_gaze.RightEye.GazeOrigin.InTrackBoxCoordinateSystem(1:2) > 0.15);
            if left_validity && right_validity
                validityColor = [0 1 0];
            end
        end
        
        origin = [screenXpixels/4 screenYpixels/4];
        size = [screenXpixels/2 screenYpixels/2];
        
        penWidthPixels = 3;
        baseRect = [0 0 size(1) size(2)];
        frame = CenterRectOnPointd(baseRect, screenXpixels/2, yCenter);
        
        Screen('FrameRect', window, validityColor, frame, penWidthPixels);
        
        % Left Eye
        if last_gaze.LeftEye.GazeOrigin.Validity.value % edited Adam
            distance = [distance; round(last_gaze.LeftEye.GazeOrigin.InUserCoordinateSystem(3)/10,1)];
            left_eye_pos_x = double(1-last_gaze.LeftEye.GazeOrigin.InTrackBoxCoordinateSystem(1))*size(1) + origin(1);
            left_eye_pos_y = double(last_gaze.LeftEye.GazeOrigin.InTrackBoxCoordinateSystem(2))*size(2) + origin(2);
            Screen('DrawDots', window, [left_eye_pos_x left_eye_pos_y], dotSizePix, validityColor, [], 2);
        end
        
        % Right Eye
        if last_gaze.RightEye.GazeOrigin.Validity.value % edited Adam
            distance = [distance;round(last_gaze.RightEye.GazeOrigin.InUserCoordinateSystem(3)/10,1)];
            right_eye_pos_x = double(1-last_gaze.RightEye.GazeOrigin.InTrackBoxCoordinateSystem(1))*size(1) + origin(1);
            right_eye_pos_y = double(last_gaze.RightEye.GazeOrigin.InTrackBoxCoordinateSystem(2))*size(2) + origin(2);
            Screen('DrawDots', window, [right_eye_pos_x right_eye_pos_y], dotSizePix, validityColor, [], 2);
        end
        pause(0.05);
    end
    
    DrawFormattedText(window, sprintf('Current distance to the eye tracker: %.2f cm.',mean(distance)), 'center', screenYpixels * 0.85, white);
    
    
    % Flip to the screen. This command basically draws all of our previous
    % commands onto the screen.
    % For help see: Screen Flip?
    Screen('Flip', window);
    
end
eyetracker.stop_gaze_data();

%% Calibrate
KbReleaseWait; % waits until all keys on the keyboard are released
while ~KbCheck % 
    Screen('TextSize', window, 20);
    DrawFormattedText(window, 'Press any key to start the calibration.', 'center', screenYpixels * 0.5, white);
    Screen('Flip', window);
end

spaceKey = KbName('Space');
RKey = KbName('R');

% Specify calibration dots on the screen
dotSizePix = 30;
dotColor = [[1 0 0];[1 1 1]]; % Red and white
leftColor = [1 0 0]; % Red
rightColor = [0 0 1]; % Bluesss

% Calibration points
lb = 0.1;  % left bound
xc = 0.5;  % horizontal center
rb = 0.9;  % right bound
ub = 0.1;  % upper bound
yc = 0.5;  % vertical center
bb = 0.9;  % bottom bound

points_to_calibrate = [[lb,ub];[rb,ub];[xc,yc];[lb,bb];[rb,bb]];

% Create calibration object
calib = ScreenBasedCalibration(eyetracker);

calibrating = true;

while calibrating
    % Enter calibration mode
    calib.enter_calibration_mode();

    for i=1:length(points_to_calibrate)

        Screen('DrawDots', window, points_to_calibrate(i,:).*screen_pixels, dotSizePix, dotColor(1,:), [], 2);
        Screen('DrawDots', window, points_to_calibrate(i,:).*screen_pixels, dotSizePix*0.5, dotColor(2,:), [], 2);

        Screen('Flip', window);

        % Wait a moment to allow the user to focus on the point
        pause(1);

        if calib.collect_data(points_to_calibrate(i,:)) ~= CalibrationStatus.Success
            % Try again if it didn't go well the first time.
            % Not all eye tracker models will fail at this point, but instead fail on ComputeAndApply.
            calib.collect_data(points_to_calibrate(i,:));
        end

    end

    DrawFormattedText(window, 'Calculating calibration result....', 'center', 'center', white);

    Screen('Flip', window);

    % Blocking call that returns the calibration result
    calibration_result = calib.compute_and_apply();

    calib.leave_calibration_mode();

    % print calibration status to command window    
    disp(' ')
    disp('Calibration Results:')
    disp('-------------------')
    if CalibrationStatus.Success==1
        disp('Calibration Successful.')
    else
        disp('Calibration Unsuccessful.')        
    end
        
    
    if calibration_result.Status ~= CalibrationStatus.Success
        break
    end
    

    
    % Calibration Result

    points = calibration_result.CalibrationPoints;

    for i=1:length(points)
        Screen('DrawDots', window, points(i).PositionOnDisplayArea.*screen_pixels, dotSizePix*0.5, dotColor(2,:), [], 2);
        for j=1:length(points(i).RightEye)
            if points(i).LeftEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                Screen('DrawDots', window, points(i).LeftEye(j).PositionOnDisplayArea.*screen_pixels, dotSizePix*0.3, leftColor, [], 2);
                Screen('DrawLines', window, ([points(i).LeftEye(j).PositionOnDisplayArea; points(i).PositionOnDisplayArea].*screen_pixels)', 2, leftColor, [0 0], 2);
            end
            if points(i).RightEye(j).Validity == CalibrationEyeValidity.ValidAndUsed
                Screen('DrawDots', window, points(i).RightEye(j).PositionOnDisplayArea.*screen_pixels, dotSizePix*0.3, rightColor, [], 2);
                Screen('DrawLines', window, ([points(i).RightEye(j).PositionOnDisplayArea; points(i).PositionOnDisplayArea].*screen_pixels)', 2, rightColor, [0 0], 2);
            end
        end

    end

    DrawFormattedText(window, 'Press the ''R'' key to recalibrate or ''Space'' to continue....', 'center', screenYpixels * 0.95, white);

    Screen('Flip', window);

  %  break;
    while 1.
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        keyCode = find(keyCode, 1);

        if keyIsDown
            if keyCode == spaceKey
                calibrating = false;
                break;
            elseif keyCode == RKey
                break;
            end
            KbReleaseWait;
        end
    end
    
end


%% Validate eyetracker calibration
KbReleaseWait; % waits until all keys on the keyboard are released
while ~KbCheck % 
    Screen('TextSize', window, 20);
    DrawFormattedText(window, 'Press any key to start the validation.', 'center', screenYpixels * 0.5, white);
    Screen('Flip', window);
end

spaceKey = KbName('Space');
RKey = KbName('R');

% Specify validation dots on the screen
dotSizePix = 30;
dotColor = [[0 1 0];[1 1 1]]; % green and white
leftColor = [1 0 0]; % Red
rightColor = [0 0 1]; % Bluesss
 
% Validation points
lb = 0.1;  % left bound
xc = 0.5;  % horizontal center
rb = 0.9;  % right bound
ub = 0.1;  % upper bound
yc = 0.5;  % vertical center
bb = 0.9;  % bottom bound 
points_to_collect = [[lb,ub];[rb,ub];[xc,yc];[lb,bb];[rb,bb]];

% Create validation object
sample_count = 30;
time_out_ms = 1000;
calib = ScreenBasedCalibrationValidation(eyetracker, sample_count, time_out_ms);

validating = true;
while validating
    calib.enter_validation_mode();
    
    for i=1:length(points_to_collect)
        calib.collect_data(points_to_collect(i,:));
        
        Screen('DrawDots', window, points_to_collect(i,:).*screen_pixels, dotSizePix, dotColor(1,:), [], 2);
        Screen('DrawDots', window, points_to_collect(i,:).*screen_pixels, dotSizePix*0.5, dotColor(2,:), [], 2);
        
        Screen('Flip', window);
        
        % Wait a moment to allow the user to focus on the point
        pause(1); 
        
    end
    
    valid_result = calib.compute();
    calib.leave_validation_mode();
    
    Screen('Flip', window);
    
   % Maybe show results on the screen from the 'valid_result' variable ? 
   % DrawFormattedText(window, 'Results.....', 'center', screenYpixels * 0.1, white);
    
   DrawFormattedText(window, 'Validation Done.', 'center', screenYpixels * 0.5, white); 
   DrawFormattedText(window, 'Press the ''R'' key to re-run validation or ''Space'' to continue....', 'center', screenYpixels * 0.95, white);
    
    Screen('Flip', window);
     

    while 1.
        [ keyIsDown, seconds, keyCode ] = KbCheck;
        keyCode = find(keyCode, 1);
        
        if keyIsDown
            if keyCode == spaceKey
                validating = false;
                break;
            elseif keyCode == RKey
                break;
            end
            KbReleaseWait;
        end
    end
    
end

%% Validation results is in "valid_result" variable:

% AverageAccuracy
% The accuracy in degrees averaged over all collected points for the left eye.

% AveragePrecision
% The precision (standard deviation) in degrees averaged over all collected points for the left eye.

 % AveragePrecisionRMS
 % The precision (root mean square of sample-to-sample error) in degrees averaged over all collected points
 % for the left eye.


disp(' ')
disp('Validation Results:')
disp('-------------------')
disp(['Average Accuracy LeftEye: ' num2str(valid_result.AverageAccuracyLeftEye)]); 
disp(['Average Precision LeftEye: ' num2str(valid_result.AveragePrecisionLeftEye)]); 
disp(['Average Precision RMS LeftEye: ' num2str(valid_result.AveragePrecisionRMSLeftEye)]); 
disp(' ');
disp(['Average Accuracy RightEye: ' num2str(valid_result.AverageAccuracyRightEye)]); 
disp(['Average Precision RightEye: ' num2str(valid_result.AveragePrecisionRightEye)]); 
disp(['Average Precision RMS RightEye: ' num2str(valid_result.AveragePrecisionRMSRightEye)]); 
 

%% Simple data colecting experiment
% pointCount = 5;
% 
% % Generate an array with coordinates of random points on the display area
% rng;
% points = rand(pointCount,2);
%
% % Start to collect data
% eyetracker.get_gaze_data();
%
% collection_time_s = 2; % seconds
%
% % Cell array to store events
% events = cell(2, pointCount);
%
% for i=1:pointCount
%     Screen('DrawDots', window, points(i,:).*screen_pixels, dotSizePix, dotColor(2,:), [], 2);
%     Screen('Flip', window); 
%     % Event when startng to show the stimulus
%     events{1,i} = {Tobii.get_system_time_stamp, points(i,:)};
%     pause(collection_time_s);
%     % Event when stopping to show the stimulus
%     events{2,i} = {Tobii.get_system_time_stamp, points(i,:)};
% end
%
% % Retreive data collected during experiment
% collected_gaze_data = eyetracker.get_gaze_data();
%
% eyetracker.stop_gaze_data();

%%  Clear
%Clear the screen. "sca" is short hand for "Screen CloseAll". This clears all features related to PTB. Note: we leave the variables in the workspace so you can have a look at them if you want. For help see: help sca
KbReleaseWait; % waits until all keys on the keyboard are released
sca;