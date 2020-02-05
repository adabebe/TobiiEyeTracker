%% ver Nov 12

clear
close all
% view info in xdf files
%% Load EEG data
% %old test
% EEG=pop_biosig(['Test data2/TestEyeTrigs.bdf']);
% test_streams = load_xdf(['Test data2/sub-P001_ses-S001_task-T1_run-001_eeg.xdf']);

% %12/11 #1
% EEG=pop_biosig(['t12_11_test/eyetracking_12_11_ver0.bdf']);
% test_streams = load_xdf(['t12_11_test/sub-P001_ses-S001_task-T1_run-001_eeg.xdf']);

%12/11 #2
EEG=pop_biosig(['t12_11_test/eyetracking_12_11_ver1.bdf']);
test_streams = load_xdf(['t12_11_test/sub-P001_ses-S001_task-T1_run-002_eeg.xdf']);

%% Parameters
sync_trigger=98; % trigger that is in the EEG signal as well in LSL/eyetracking data - trial onset
sync_trigger_off=99; % trigger that is in the EEG signal as well in LSL/eyetracking data - trial offset
eog_chan=131; % EEG/EOG channel corresponding to horizontal eye movements 
left_mastoid_chan=129;
right_mastoid_chan=130;

%% re-reference EOG/EEG to mastoids
% not necessary but leads to better data
if (~isempty(left_mastoid_chan))&&(~isempty(right_mastoid_chan))
    mastoids_avg=mean(EEG.data([left_mastoid_chan right_mastoid_chan],:));
    tmp = EEG.data';
    refEEG = tmp-repmat(mastoids_avg',[1,size(EEG.data,1)]); 
    EEG.data=refEEG';    
end

%% Parse LSL/Tobii data
tobii_stream_idx=[]; pres_stream_idx=[]; eeg_stream_idx=[];
% access name of each stream
for i = 1:size(test_streams,2)
    stream_name = test_streams{1, i}.info.name;
    if strcmp(stream_name,'Tobii')
        tobii_stream_idx = i;
    elseif strcmp(stream_name,'Presentation')
        pres_stream_idx = i;
    elseif strcmp(stream_name,'Biosemi')
        eeg_stream_idx = i;
    else
        warning('Unknown stream found in xdf')
    end
end
if ~isempty(tobii_stream_idx) % Get Tobii data
    tobii_data = test_streams{1, tobii_stream_idx}.time_series;
    tobii_eye_coords = tobii_data([1:2,4:5],:);
    tobii_pupil_data = tobii_data([3,6],:);
    tobii_time  = test_streams{1, tobii_stream_idx}.time_stamps;
end
lsl_trig_code=[];
lsl_trig_time=[];
if ~isempty(pres_stream_idx) % Get presentation data
    pres_cell_events = test_streams{1, pres_stream_idx}.time_series;
    for i=1:size(pres_cell_events,2)
        % keep only event numbers/codes
        curr_event = pres_cell_events{1,i};
        tf = strfind(curr_event,{'ecode'});
        isolated_event_code = curr_event(tf(1):tf(2));
        % keep numbers only
        curr_e_string = isolated_event_code(and(double(isolated_event_code)>47,double(isolated_event_code)<58));
        if ~isempty(curr_e_string) % ignore empty 'quit' trigger as it does not show in .bdf
            lsl_trig_code(i) = str2num(curr_e_string);
            lsl_trig_time(i)  = test_streams{1, pres_stream_idx}.time_stamps(i);
        end
    end
end

%% plot LSL data (tobii) along with triggers
% figure(3)
% hold on
% plot(tobii_time-tobii_time(1), tobii_eye_coords(1,:),'k'),
% plot(tobii_time-tobii_time(1), tobii_eye_coords(2,:),'k'),
% plot(tobii_time-tobii_time(1), tobii_eye_coords(3,:),'k'),
% plot(tobii_time-tobii_time(1), tobii_eye_coords(4,:),'k'),
% for i=1:length(lsl_trig_time)
%     line([lsl_trig_time(i)-tobii_time(1),lsl_trig_time(i)-tobii_time(1)],[0,1.2],'color','g')
%     text(lsl_trig_time(i),1,num2str(lsl_trig_code(i)),'HorizontalAlignment','left','color','g');
% end
% hold off

%% Fill missing Tobii eye values
% TODO maybe keep track of interpolated values ??
for i=1:4
    tobii_eye_coords(i,:) = fillmissing(tobii_eye_coords(i,:),'linear');
end

%% Resample Tobii to the rate of EEG.srate
tobii_time_rs=min(tobii_time):(1/EEG.srate):max(tobii_time);
tobii_eye_coords_rs=[];
for i=1:size(tobii_eye_coords,1)
    tobii_eye_coords_rs(i,:)=interp1(tobii_time,tobii_eye_coords(i,:),tobii_time_rs);
end


%% shift LSL data to be in sync with EEG based using ONSET trigger only
tobii_time_matched=[];tobii_eye_coords_matched=[];lsl_trig_time_matched=[]; lsl_trig_code_matched=[];
sync_trigger_LSL_indxs=find(lsl_trig_code==sync_trigger);
tobii_time_rs_old=tobii_time_rs;
lsl_trig_time_old=lsl_trig_time;
try
    for j=1:length(sync_trigger_LSL_indxs) % for each sync trigger in LSL
        tobii_time_rs=tobii_time_rs_old;
        lsl_trig_time=lsl_trig_time_old;
        
        %% find where is the corresponding sync trigger in EEG
        sync_trig_eeg_indxs=find([EEG.event(:).type]==sync_trigger);
        if length(sync_trig_eeg_indxs)~=length(sync_trigger_LSL_indxs); error('WARNING: Different number of sync triggers in LSL stream and EEG file !!'); end
        sync_trig_eeg_latency_s=(EEG.event(sync_trig_eeg_indxs(j)).latency-1)./EEG.srate; % dont forget -1 because the first sample point corresponds to time 0; see https://sccn.ucsd.edu/wiki/Chapter_03:_Event_Processing
        
        %% find where is the sync trigger in tobii and shift the time accordingly so it matches EEG
        [~ ,idx] = min(abs(tobii_time_rs - lsl_trig_time(sync_trigger_LSL_indxs(j))));
        tobii_time_rs=tobii_time_rs-tobii_time_rs(idx)+sync_trig_eeg_latency_s;  % shift the Tobii time so that the sync trigger has the time that corresponds to the sync trigger time in EEG
        lsl_trig_time_rs=lsl_trig_time-lsl_trig_time(sync_trigger_LSL_indxs(j))+sync_trig_eeg_latency_s; % shift the LSL trigger time correspondingly...
        
        %% cut out the Tobii signal where EEG was not saving (not between 98-99 trigger)
        sync_trig_off_eeg_indxs=find([EEG.event(:).type]==sync_trigger_off);   % find all sync offset triggers in EEG
        sync_trig_off_eeg_indxs=sync_trig_off_eeg_indxs(sync_trig_off_eeg_indxs>sync_trig_eeg_indxs(j)); %find offset trigger that is after current sync onset trigger; => sync_trig_off_eeg_indxs(1) is the one
        n_samples_betw_sync_trigs_eeg=round(EEG.event(sync_trig_off_eeg_indxs(1)).latency-EEG.event(sync_trig_eeg_indxs(j)).latency); % distance in samples between 98 and 99 trigs in EEG
        % as EEG event samples can be a fraction (if downsampled, it is better to run this before EEG downsampling), this should to be rounded        
        tobii_time_matched=[tobii_time_matched tobii_time_rs(idx:(idx+n_samples_betw_sync_trigs_eeg-1))];
        tobii_eye_coords_matched=[tobii_eye_coords_matched tobii_eye_coords_rs(:,idx:idx+n_samples_betw_sync_trigs_eeg-1)];
        
        
        %% cut out LSL triggers where EEG was not saving (not between 98-99 trigger)
        sync_trig_eeg_latency_s_off=(EEG.event(sync_trig_off_eeg_indxs(1)).latency-1)./EEG.srate; % dont forget -1 because the first sample point corresponds to time 0; see https://sccn.ucsd.edu/wiki/Chapter_03:_Event_Processing
        lsl_trig_times_to_keep=lsl_trig_time_rs(and((lsl_trig_time_rs>=sync_trig_eeg_latency_s),(lsl_trig_time_rs<=ceil(sync_trig_eeg_latency_s_off))));
        lsl_trig_codes_to_keep=lsl_trig_code(and((lsl_trig_time_rs>=sync_trig_eeg_latency_s),(lsl_trig_time_rs<=ceil(sync_trig_eeg_latency_s_off))));        
        lsl_trig_time_matched=[lsl_trig_time_matched lsl_trig_times_to_keep];
        lsl_trig_code_matched=[lsl_trig_code_matched lsl_trig_codes_to_keep];
        
    end
    
    %% zero-pad: there are EEG data samples before and after sync trigger, so zero pad the Tobii samples => tobii and EEG share time axes from now on
    EEG_times_in_sec=EEG.times/1000;
    tobii_eye_coords_matched2=nan(size(tobii_eye_coords_matched,1),EEG.pnts);
    tobii_eye_coords_matched2(:,ismember(EEG_times_in_sec,tobii_time_matched))=tobii_eye_coords_matched(:,:);
    tobii_time_matched2=1e-3*EEG.times;
catch
    error('Cannot sync data. Are all sync triggers present in EEG and Tobii?');
end

%% For debugging purposes only !- introduce artficial lag into eye-tracking => should be reflected in results ! .. should be commented out if not debugging
% test_shift=3;
% for i=1:size(tobii_eye_coords_matched2,1)
%     tobii_eye_coords_matched2(i,:)=[tobii_eye_coords_matched2(i,test_shift:end) zeros(1,test_shift-1)];
% end

%% Plot EOG and TObii
figure(1)
hold on

% Tobii Signals 
plot(tobii_time_matched2, tobii_eye_coords_matched2(1,:),'r'),  % plot only left horizontal
%plot(tobii_time_matched2, tobii_eye_coords_matched2(2,:),'k'),
%plot(tobii_time_matched2, tobii_eye_coords_matched2(3,:),'k'),
%plot(tobii_time_matched2, tobii_eye_coords_matched2(4,:),'k'),

for i=1:length(lsl_trig_code_matched) % Presentation Triggers
    line([lsl_trig_time_matched(i),lsl_trig_time_matched(i)],[0,1.2],'color','g')
    text(lsl_trig_time_matched(i),1,num2str(lsl_trig_code_matched(i)),'HorizontalAlignment','left','color','g');
end

for i=1:size(EEG.event,2)  % EEG Triggers
    line([1e-3*EEG.times(round(EEG.event(i).latency)),1e-3*EEG.times(round(EEG.event(i).latency))],[0, -1.2])
    text(1e-3*EEG.times(round(EEG.event(i).latency)),-1,num2str(EEG.event(i).type),'HorizontalAlignment','left','color','b');
end

EOG=EEG.data(eog_chan,:); % Biosemi EOG/EEG
EOG=EOG-mean(EOG);
EOG=EOG/(max(EOG)-min(EOG));
plot(1e-3*EEG.times, -EOG)
hold off

%% calculate cross correlation between eye muscle data and eye tracking data

MAXLAG=EEG.srate*0.5; % run xcorr for 0.5 sec lags
sync_trig_eeg_indxs=find([EEG.event(:).type]==sync_trigger);
sync_trig_off_eeg_indxs=find([EEG.event(:).type]==sync_trigger_off);
trial_length_sec=floor(((EEG.event(sync_trig_off_eeg_indxs(1)).latency-1)./EEG.srate)-((EEG.event(sync_trig_eeg_indxs(1)).latency-1)./EEG.srate)); % just for corr calculation

EOG=EEG.data(eog_chan,:); % Scaled Biosemi EOG/EEG
EOG=EOG-mean(EOG);
EOG=-EOG/(max(EOG)-min(EOG));

% Correlate only epoched data without transitions
sync_trig_eeg_indxs=find([EEG.event(:).type]==sync_trigger);   % find all sync offset triggers in EEG
xc=[];
for trial=1:size(sync_trig_eeg_indxs,2)
    trial_indxs=(EEG.event(sync_trig_eeg_indxs(trial)).latency):(EEG.event(sync_trig_eeg_indxs(trial)).latency)+trial_length_sec*EEG.srate;
    EOGforCorr=EOG(:,trial_indxs);
    TobiiforCorr=tobii_eye_coords_matched2(1,trial_indxs);
    [xc(trial,:),lags] = xcorr(zscore(EOGforCorr),zscore(TobiiforCorr),MAXLAG);
  %  plot(zscore(EOGforCorr)); hold on; plot(zscore(TobiiforCorr));
end
mean_xc=mean(xc);
[~,indx]=max(abs(mean_xc));
time_lag_mesured=lags(indx)./EEG.srate*1000;

%% Plot results
figure(2)
subplot(1,2,1)
time_lag_ms=1000*lags./EEG.srate;
plot(time_lag_ms,mean_xc)
hold on
ylims=get(gca,'ylim');
plot([time_lag_mesured time_lag_mesured],[ylims(1) max(mean_xc)])

xlabel('time lag (ms)')
ylabel('cross-correlation')
title({'Cross-correlation EOG vs EYE tracker',['Lag= ' num2str(time_lag_mesured) ' (ms)']})

subplot(1,2,2)
plot(time_lag_ms,mean_xc)
xlim([-25 25])
hold on
ylims=get(gca,'ylim');
plot([time_lag_mesured time_lag_mesured],ylims)
xlabel('time lag (ms)')
ylabel('cross-correlation')
title({'Cross-correlation EOG vs EYE tracker (zoomed)',['Lag= ' num2str(time_lag_mesured) ' (ms)']})

disp('**** Eyetracking co-registration results: *****')
disp(' ')
disp(['Cross-correlation lag between EOG and Eyetracking is ' num2str(time_lag_mesured) ' (ms)/ ' num2str(lags(indx)) ' (samples)'])
%% Calculate Timing error based on offset sync triggers in Tobii/LSL and EEG
sync_trig_off_eeg_indxs=find([EEG.event(:).type]==sync_trigger_off);   % find all sync offset triggers in EEG
sync_trig_off_eeg_latency_s=([EEG.event(sync_trig_off_eeg_indxs).latency]-1)./EEG.srate;
sync_trig_off_lsl_indxs=find(lsl_trig_code_matched==sync_trigger_off);
sync_trig_off_lsl_latency_s=lsl_trig_time_matched(sync_trig_off_lsl_indxs);

mean_diff=mean(abs((sync_trig_off_eeg_latency_s-sync_trig_off_lsl_latency_s)))*1000; % ms
max_diff=max(abs((sync_trig_off_eeg_latency_s-sync_trig_off_lsl_latency_s)))*1000; % ms
min_diff=min(abs((sync_trig_off_eeg_latency_s-sync_trig_off_lsl_latency_s)))*1000; % ms
std_diff=std(abs((sync_trig_off_eeg_latency_s-sync_trig_off_lsl_latency_s)))*1000; % ms
disp(' ')
disp(['Timing Error based on offset triggers (n=' num2str(length(sync_trig_off_lsl_indxs)) '):'])
disp(['Mean: ' num2str(mean_diff) '(ms)'])
disp(['Max: ' num2str(max_diff) '(ms)'])
disp(['Min: ' num2str(min_diff) '(ms)'])
disp(['Std: ' num2str(std_diff) '(ms)'])
