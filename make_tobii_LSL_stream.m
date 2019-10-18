%% These examples show how to transmit a numeric multi-channel time series through LSL

addpath(genpath('../liblsl-Matlab'))

%% Sample-by-Sample make test stream;
% https://github.com/labstreaminglayer/liblsl-Matlab/blob/master/examples/SendData.m

% instantiate the library
disp('Loading library...');
lib = lsl_loadlib();

% make a new stream outlet
disp('Creating a new streaminfo...');
info = lsl_streaminfo(lib,'Tobii','EyeTracker',8,100,'cf_float32','sdfwerr32432');

disp('Opening an outlet...');
outlet = lsl_outlet(info);

% send data into the outlet, sample by sample
disp('Now transmitting data...');
while true
    outlet.push_sample(randn(8,1));
    pause(0.01);
end



