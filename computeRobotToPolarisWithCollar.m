% Compute robot to polaris transform by comparing the endoscope tip
% position as reported by the API (kinematics) against the observed motion
% of the xi tracking collar
%
% Can be run:
% * At specified observations indicated in captureNoteList
% * Between specified start and stop times w.r.t. synchronization start time
% * At specified observations subject to start and stop times

function [rmse,TF_robot_to_polaris] = computeRobotToPolarisWithCollar(...
    polarisFile, ...
    kinematicsFile, ...
    collarToolID, ...
    scopeArmID, ...
    varargin)

% accept additional arguments:
useStartStopTimes = false;
useCaptureNoteList = false;
switch(nargin)
    case 5  % ... captureNoteList);
        captureNoteList = varargin{1};
        useCaptureNoteList = true;
    case 6  % ... startTimeSec, stopTimeSec);
        startTimeSec = varargin{1};
        stopTimeSec = varargin{2};
        useStartStopTimes = true;
    case 7 % ... captureNoteList, startTimeSec, stopTimeSec);
        captureNoteList = varargin{1};
        startTimeSec = varargin{2};
        stopTimeSec = varargin{3};
        useCaptureNoteList = true;
        useStartStopTimes = true;
end

% load kinematics file
jv_tab = readtable(kinematicsFile,'NumHeaderLines',0,'Delimiter',',');
jv_ecm_mask = strcmp(jv_tab.Var3,['ISI_USM' num2str(scopeArmID)]);  % NOTE: NEED TO DEFINE WHICH ARM SCOPE IS IN!
jv_ecm_subtab = jv_tab(jv_ecm_mask,:);
allServoTimes = table2array(jv_ecm_subtab(:,2));
allServoTimes = unwrap( allServoTimes*2*pi/65535 )*65535/(2*pi);
allServoTimes = allServoTimes - allServoTimes(1);
ecm_sec = allServoTimes*(7.5e-4);     % DISCOVERED THIS EMPIRICALLY (750us/count = 1333.33Hz), may want to check

% load polaris file and adjust polaris time
polaris_tab = readtable(polarisFile);
polaris_collar_mask = strcmp(polaris_tab.tool_id,collarToolID);
subtab = polaris_tab(polaris_collar_mask,:);
polaris_sec = subtab.polaris_time - (polaris_tab(1,:).polaris_time );

% if using a capture list restrict the polaris subtab to just those
% specific observations
if(useCaptureNoteList)
    captureNoteMask = zeros(size(subtab,1),1);
    for pose_idx = 1:length(captureNoteList)
        thisCaptureNote = captureNoteList{pose_idx};
        thisNoteMask = strcmp(thisCaptureNote,subtab.capture_note);
        switch(nnz(thisNoteMask))
            case 0
                warning('No observations of ''%s''...',thisCaptureNote);
            case 1
                captureNoteMask = captureNoteMask | thisNoteMask;
            otherwise
                warning('More than one observation of ''%s''...',thisCaptureNote);
        end
    end
    subtab(~captureNoteMask,:) = [];
    polaris_sec(~captureNoteMask,:) = [];
end

% for each collar observation find the closest ECM kinematic datapoint
cal_kinpts = nan(3,size(subtab,1));
cal_TF = nan(4,4,size(subtab,1));
actual_ecm_sec = nan(size(polaris_sec));
for i = 1:size(subtab,1)
    [minval,ecm_idx] = min(abs(ecm_sec - polaris_sec(i)));
    if(minval < 0.100)
        actual_ecm_sec(i) = ecm_sec(ecm_idx);
        TFjv = deserialize_isi_tf(table2array(jv_ecm_subtab(ecm_idx,4:15)));
        TFjv(1:3,4) = TFjv(1:3,4)*1000; % convert to millimeters
        cal_kinpts(:,i) = TFjv(1:3,4);
        q = [subtab(i,:).q0 subtab(i,:).q1 subtab(i,:).q2 subtab(i,:).q3];
        t = [subtab(i,:).tx subtab(i,:).ty subtab(i,:).tz];
        cal_TF(:,:,i) = [quat2matrix(q), t'; 0 0 0 1];
    else
        error('No suitable match for Polaris observation #%d, min error: %0.3fsec',i,minval);
    end
end
nanmask = isnan(cal_kinpts(1,:));
if(useStartStopTimes)
    nanmask = nanmask & ~((polaris_sec > startTimeSec*60) & (polaris_sec < stopTimeSec*60))';
end
cal_kinpts(:,nanmask) = [];
cal_TF(:,:,nanmask) = [];


% run constrained optimization, restricting search space to 
x0 = [0, 0, -600]';  % initial guess for [x,y,z]' coordinates of scope tip in collar frame
f = @(x_current)robot_to_polaris_rmse(x_current,cal_kinpts,cal_TF);
options = optimset('MaxFunEvals',1e8);
A = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
b = [10 10 10 10 -400 800];
[x_opt, ~] = fmincon(f,x0,A,b,[],[],[],[],[],options);
[finalRMSE,finalTF] = robot_to_polaris_rmse(x_opt,cal_kinpts,cal_TF);
rmse = finalRMSE;
TF_robot_to_polaris = finalTF;

% display results for debugging
figure; 
hold on; grid on; axis equal;
cal_kinpts_polaris_frame = hTF(cal_kinpts,finalTF,0);
ph_kinematics = plot3(cal_kinpts_polaris_frame(1,:),cal_kinpts_polaris_frame(2,:),cal_kinpts_polaris_frame(3,:),'o','MarkerSize',10,'LineWidth',2,'Color',[0.8 0 0.8]);
cal_colpts = [];
for i = 1:size(cal_TF,3)
    tippt = hTF(x_opt,cal_TF(:,:,i),0);
    %     plot3(tippt(1),tippt(2),tippt(3),'-','LineWidth',2,'MarkerSize',5,'Color',[0.8 0 0]);
    cal_colpts(:,end+1) = tippt;
end
ph_polaris = plot3(cal_colpts(1,:),cal_colpts(2,:),cal_colpts(3,:),'.','MarkerSize',30,'Color',[0 0 0.8]);
view([-81,16.4]);
xlim([-320,-50]);
ylim([-120,80]);
zlim([-2280,-2080]);
legend([ph_polaris,ph_kinematics],{'Polaris','Kinematics'});
end


% objective function for finding T_robot_to_polaris using tracking collar
function [rmse,TF_robot_to_polaris] = robot_to_polaris_rmse(x_current,cal_kinpts,cal_TF)
cal_colpts = nan(size(cal_kinpts));
for data_idx = 1:size(cal_TF,3)
    TF_collar_to_polaris = cal_TF(:,:,data_idx);
    cal_colpts(:,data_idx) = hTF(x_current,TF_collar_to_polaris,0);
end
[newKP,TF_robot_to_polaris,RMSE] = rigid_align_svd(cal_kinpts,cal_colpts);
rmse = RMSE;
end
