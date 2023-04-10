% QTM Analysis for One Subject 
% Clear workspace and command window
% clear,clc, close all 

SubjectID = 'PSCI06 SLA and PSCI07 SB';

%% Preallocate cell arrays
QTMData.SLA = cell(1,6);
QTMData.SB = cell(1,6);
QTMFilt.SLA = cell(1,6);
QTMFilt.SB = cell(1,6);

%% Pull load cell and marker trajectory data into a struct

% Create cell array of the marker trajectory data 
TempQTM.SLA = {SLA_0001,SLA_0002,SLA_0003,SLA_0004,SLA_Baseline,SLA_Retention};
TempQTM.SB = {SB_0001,SB_0002,SB_0003,SB_0004,SB_Baseline,SB_Retention};

for i=1:6 
    
    % Pull each channel of the load cell data and add it to the QTMData struct.
    % SLA files
    QTMData.SLA{1,i}.Load = array2table(squeeze(TempQTM.SLA{1, i}.Analog.Data'));           
    % SB files 
    QTMData.SB{1,i}.Load = array2table(squeeze(TempQTM.SB{1, i}.Analog.Data'));

    % Sum channels 1:4 to get Right load cell data. Sum channels 5:8 to get Left load cell data
    QTMData.SLA{1,i}.Load.RLoad = sum(QTMData.SLA{1,i}.Load{:,1:4},2);
    QTMData.SLA{1,i}.Load.LLoad = sum(QTMData.SLA{1,i}.Load{:,5:8},2); 
    QTMData.SB{1,i}.Load.RLoad = sum(QTMData.SB{1,i}.Load{:,1:4},2);
    QTMData.SB{1,i}.Load.LLoad = sum(QTMData.SB{1,i}.Load{:,5:8},2);

    % Pull each dimenstion (X,Y,Z) of the trajectory data and make a struct 
    % SLA files
        QTMData.SLA{1,i}.Trajectories = struct('X',array2table(squeeze(TempQTM.SLA{1, i}.Trajectories.Labeled.Data(:,1,:))',...
        'VariableNames',{'R_IC','L_IC','R_GTO','L_GTO','R_KNEE','L_KNEE','R_ANK','L_ANK','R_5MTP','L_5MTP'}),...
        'Y',array2table(squeeze(TempQTM.SLA{1, i}.Trajectories.Labeled.Data(:,2,:))',...
        'VariableNames',{'R_IC','L_IC','R_GTO','L_GTO','R_KNEE','L_KNEE','R_ANK','L_ANK','R_5MTP','L_5MTP'}),...
        'Z',array2table(squeeze(TempQTM.SLA{1, i}.Trajectories.Labeled.Data(:,3,:))',...
        'VariableNames',{'R_IC','L_IC','R_GTO','L_GTO','R_KNEE','L_KNEE','R_ANK','L_ANK','R_5MTP','L_5MTP'}));
    % SB files 
        QTMData.SB{1,i}.Trajectories = struct('X',array2table(squeeze(TempQTM.SB{1, i}.Trajectories.Labeled.Data(:,1,:))',...
        'VariableNames',{'R_IC','L_IC','R_GTO','L_GTO','R_KNEE','L_KNEE','R_ANK','L_ANK','R_5MTP','L_5MTP'}),...
        'Y',array2table(squeeze(TempQTM.SB{1, i}.Trajectories.Labeled.Data(:,2,:))',...
        'VariableNames',{'R_IC','L_IC','R_GTO','L_GTO','R_KNEE','L_KNEE','R_ANK','L_ANK','R_5MTP','L_5MTP'}),...
        'Z',array2table(squeeze(TempQTM.SB{1, i}.Trajectories.Labeled.Data(:,3,:))',...
        'VariableNames',{'R_IC','L_IC','R_GTO','L_GTO','R_KNEE','L_KNEE','R_ANK','L_ANK','R_5MTP','L_5MTP'}));
end 


% Turn data into table, adding trial labels
QTMData.SLA = cell2table(QTMData.SLA, 'VariableNames',{'SLA_1', 'SLA_2', 'SLA_3', 'SLA_4',...
    'SLA_Baseline', 'SLA_Retention'});

QTMData.SB = cell2table(QTMData.SB, 'VariableNames',{'SB_1', 'SB_2', 'SB_3', 'SB_4',...
    'SB_Baseline', 'SB_Retention'});


%% Filtering (Lowpass_Test function in Lab_Matlab_Scripts) 
% inputs are a vector of QTM data, the sampling rate, and a cutoff
% QTMFilt = Lowpass_Test(QTMData.SLA.Load, 120, CUTOFF);
% The Lowpass_Test function is different from the lowpass function because it can handle NaN value
cd('/Users/morganlea/Documents/MATLAB/Lab_Matlab_Scripts')

for i=1:6 
    % Filter Load cell data 
    % 1000 sampling rate, 20Hz cutoff for Load Cell Data
    QTMFilt.SLA{1,i}.Load = Lowpass_Test(table2array(QTMData.SLA{1,i}.Load), 1000, 20); 
    QTMFilt.SB{1,i}.Load = Lowpass_Test(table2array(QTMData.SB{1,i}.Load), 1000, 20);

    % Filter Kinematic data
    % 120 sampling rate, 6Hz cutoff for Kinematic Data (from Winters
    % Biomechanics text)
    QTMFilt.SLA{1,i}.Trajectories.X = Lowpass_Test(table2array(QTMData.SLA{1,i}.Trajectories.X), 120, 6);
    QTMFilt.SLA{1,i}.Trajectories.Y = Lowpass_Test(table2array(QTMData.SLA{1,i}.Trajectories.Y), 120, 6);
    QTMFilt.SLA{1,i}.Trajectories.Z = Lowpass_Test(table2array(QTMData.SLA{1,i}.Trajectories.Z), 120, 6);
    QTMFilt.SB{1,i}.Trajectories.X = Lowpass_Test(table2array(QTMData.SB{1,i}.Trajectories.X), 120, 6);
    QTMFilt.SB{1,i}.Trajectories.Y = Lowpass_Test(table2array(QTMData.SB{1,i}.Trajectories.Y), 120, 6);
    QTMFilt.SB{1,i}.Trajectories.Z = Lowpass_Test(table2array(QTMData.SB{1,i}.Trajectories.Z), 120, 6);

    % Add Variable Names to the Trajectory Data - Not working
    % QTMFilt.SLA{1,i}.Trajectories.X = renamevars(QTMFilt.SLA{1,i}.Trajectories.X,...
    %     ["Var1","Var2","Var3","Var4","Var5","Var6","Var7","Var8","Var9","Var10"],...
    %     ["L_IC","R_IC","L_GTO","R_GTO","L_KNEE","R_KNEE","R_ANK","L_ANK","L_5MTP","R_5MTP"]);
end 

% Turn filtered data into a table, adding trial labels. Add variable names
% to the Load and Trajectory Data 

QTMFilt.SLA = cell2table(QTMFilt.SLA, 'VariableNames',{'SLA_1', 'SLA_2', 'SLA_3', 'SLA_4',...
    'SLA_Baseline', 'SLA_Retention'});

QTMFilt.SB = cell2table(QTMFilt.SB, 'VariableNames',{'SB_1', 'SB_2', 'SB_3', 'SB_4',...
    'SB_Baseline', 'SB_Retention'});


%% Heel Strike Detection 
GaitEvents = cell(1,6);

for i = 1:6
% Define time vector to match kinetic data (i.e. 1,000Hz sampling rate)
TIME = linspace(0, height(QTMFilt.SLA{1,i}.Load)/1000., height(QTMFilt.SLA{1,i}.Load))';

% Define time vector to match kinematic data (i.e. 120 Hz sampling rate)
TIME_SPLINE = linspace(0, height(QTMFilt.SLA{1,i}.Trajectories.X)/120, height(QTMFilt.SLA{1,i}.Trajectories.X))';

% TemporalParam function
% [EVENTS, EVENTS_KINEMATICS, CADENCE, STEP_TIME, L_DST, R_DST, SW_TIME, ST_TIME]
% = TemporalParam(3, TIME, TIME_SPLINE, L_GRF_Z_FILTERED, R_GRF_Z_FILTERED, L_ANKLE, 
% R_ANKLE, FRAME_RATE, TIME_THRESH)

% Gait Events for SLA
    % Load: the 9th column is RLoad and 10th column is LLoad
    % Trajectories: the 7th column is R_Ankle and 8th column is L_Ankle
    [GaitEvents{1,i}.SLA.Kinetics,GaitEvents{1,i}.SLA.Kinematics, ~, ~, ~, ~, ~, ~] = TemporalParam...
    (3, TIME, TIME_SPLINE,  QTMFilt.SLA{1,i}.Load(:,10),  QTMFilt.SLA{1,i}.Load(:,9), ...
    QTMFilt.SLA{1,i}.Trajectories.X(:,8), QTMFilt.SLA{1,i}.Trajectories.X(:,7), 100, 50,0.25);

end 



% Gait Events for SB 

for i = 1:6
% Define time vector to match kinetic data (i.e. 1,000Hz sampling rate)
TIME = linspace(0, height(QTMFilt.SB{1,i}.Load)/1000., height(QTMFilt.SB{1,i}.Load))';

% Define time vector to match kinematic data (i.e. 120 Hz sampling rate)
TIME_SPLINE = linspace(0, height(QTMFilt.SB{1,i}.Trajectories.X)/120, height(QTMFilt.SB{1,i}.Trajectories.X))';

% TemporalParam function
% [EVENTS, EVENTS_KINEMATICS, CADENCE, STEP_TIME, L_DST, R_DST, SW_TIME, ST_TIME]
% = TemporalParam(3, TIME, TIME_SPLINE, L_GRF_Z_FILTERED, R_GRF_Z_FILTERED, L_ANKLE, 
% R_ANKLE, FRAME_RATE, TIME_THRESH)


    % Load: the 9th column is RLoad and 10th column is LLoad
    % Trajectories: the 7th column is R_Ankle and 8th column is L_Ankle
    [GaitEvents{1,i}.SB.Kinetics,GaitEvents{1,i}.SB.Kinematics, ~, ~, ~, ~, ~, ~] = TemporalParam...
    (3, TIME, TIME_SPLINE,  QTMFilt.SB{1,i}.Load(:,10),  QTMFilt.SB{1,i}.Load(:,9), ...
    QTMFilt.SB{1,i}.Trajectories.X(:,8), QTMFilt.SB{1,i}.Trajectories.X(:,7), 100, 50,0.25);

end 

% Convert gait events into a table 

for i = 1:6
GaitEvents{1,i}.SLA.Kinetics = array2table(GaitEvents{1,i}.SLA.Kinetics, ...
    'VariableNames', {'L_HS', 'L_TO', 'R_HS', 'R_TO'});
GaitEvents{1,i}.SLA.Kinematics = array2table(GaitEvents{1,i}.SLA.Kinematics, ...
    'VariableNames', {'L_HS', 'L_TO', 'R_HS', 'R_TO'});
GaitEvents{1,i}.SB.Kinetics = array2table(GaitEvents{1,i}.SB.Kinetics, ...
    'VariableNames', {'L_HS', 'L_TO', 'R_HS', 'R_TO'});
GaitEvents{1,i}.SB.Kinematics = array2table(GaitEvents{1,i}.SB.Kinematics, ...
    'VariableNames', {'L_HS', 'L_TO', 'R_HS', 'R_TO'});
end 


%% Calculate step length asymmetry

% SLA defintion by Padmanabhan et al., Journal of NeuroEngineering and
% Rehabilition, 2020
% https://link.springer.com/article/10.1186/s12984-020-00732-z

% Also from: https://journals.sagepub.com/doi/pdf/10.1177/1545968319855028
% The longer and shorter step lengths were defined
% based on the average step length of each limb over the four
% minutes of walking at self-selected speed.
% SLA = (longer step - shorter step) / (longer step + shorter step) 

% Find step length at heel strike:
% Left: 
%   L ankle at L heel strike - R ankle at L heelstrike
%   L_ANK(L_HS) - R_ANK(L_HS) 
% Right: 
%   R ankle at R heel strike - L ankle at R heelstrike
%   R_ANK(R_HS) - L_ANK(R_HS)

StepLength = cell(1,6);

for i = 1:6

% SLA
% Trajectories: the 7th column is R_Ankle and 8th column is L_Ankle 
StepLength{1,i}.SLA.Left = abs(QTMFilt.SLA{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SLA.Kinematics.L_HS,8) - QTMFilt.SLA{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SLA.Kinematics.L_HS,7))./1000;
StepLength{1,i}.SLA.Right = abs(QTMFilt.SLA{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SLA.Kinematics.R_HS,7) - QTMFilt.SLA{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SLA.Kinematics.R_HS,8))./1000;            
 

% SB
% Trajectories: the 7th column is R_Ankle and 8th column is L_Ankle            
StepLength{1,i}.SB.Left = abs(QTMFilt.SB{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SB.Kinematics.L_HS,8) - QTMFilt.SB{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SB.Kinematics.L_HS,7))./1000;
StepLength{1,i}.SB.Right = abs(QTMFilt.SB{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SB.Kinematics.R_HS,7) - QTMFilt.SB{1,i}.Trajectories.X...
    (GaitEvents{1,i}.SB.Kinematics.R_HS,8))./1000;            
 
end 


SLA_Calculated = cell(1, 6);
SLA_Matched_Strides = cell(1, 6);

for i = 1:6
% Baseline data is in column 5. 
    if mean(StepLength{1,5}.SLA.Left) > mean(StepLength{1,5}.SLA.Right)
    SLA_Calculated{1,i}.SLA = (StepLength{1,i}.SLA.Left - StepLength{1,i}.SLA.Right)./...
     (StepLength{1,i}.SLA.Left + StepLength{1,i}.SLA.Right);

    else % if baseline right step length is longer
      SLA_Calculated{1,i}.SLA = (StepLength{1,i}.SLA.Right - StepLength{1,i}.SLA.Left)./...
          (StepLength{1,i}.SLA.Right + StepLength{1,i}.SLA.Left);
    end
% For split belt, the shorter step goes at the faster speed. Faster -
% slower / (faster + slower). Roemmich 2016 "Seeing the Errors You Feel"

     if mean(StepLength{1,5}.SB.Left) > mean(StepLength{1,5}.SB.Right)
    SLA_Calculated{1,i}.SB = (StepLength{1,i}.SB.Right - StepLength{1,i}.SB.Left)./...
     (StepLength{1,i}.SB.Left + StepLength{1,i}.SB.Right);

    else % if baseline right step length is longer
      SLA_Calculated{1,i}.SB = (StepLength{1,i}.SB.Left - StepLength{1,i}.SB.Right)./...
          (StepLength{1,i}.SB.Right + StepLength{1,i}.SB.Left);
     end 
end


SLA_Calculated = cell2table(SLA_Calculated, 'VariableNames',...
    {'Trial_1','Trial_2', 'Trial_3', 'Trial_4','Baseline', 'Retention'});


SLA_timeseries = [SLA_Calculated.Baseline.SLA;SLA_Calculated.Trial_1.SLA;SLA_Calculated.Trial_2.SLA;...
    SLA_Calculated.Trial_3.SLA;SLA_Calculated.Trial_4.SLA;SLA_Calculated.Retention.SLA];


SB_timeseries = [SLA_Calculated.Baseline.SB;SLA_Calculated.Trial_1.SB;SLA_Calculated.Trial_2.SB;...
    SLA_Calculated.Trial_3.SB;SLA_Calculated.Trial_4.SB;SLA_Calculated.Retention.SB];

%% Plots

% Plot Time Series Data
figure
tiledlayout(2,1)
sgtitle(SubjectID)

nexttile
plot(SLA_timeseries);
hold on
xline(height(SLA_Calculated.Baseline.SLA),"--")
xline(height(SLA_Calculated.Baseline.SLA) + height(SLA_Calculated.Trial_1.SLA),"--")
xline(height(SLA_Calculated.Baseline.SLA) + height(SLA_Calculated.Trial_1.SLA) + height(SLA_Calculated.Trial_2.SLA),"--")
xline(height(SLA_Calculated.Baseline.SLA) + height(SLA_Calculated.Trial_1.SLA) + height(SLA_Calculated.Trial_2.SLA) ...
    + height(SLA_Calculated.Trial_3.SLA),"--")
xline(height(SLA_Calculated.Baseline.SLA) + height(SLA_Calculated.Trial_1.SLA) + height(SLA_Calculated.Trial_2.SLA) ...
    + height(SLA_Calculated.Trial_3.SLA) + height(SLA_Calculated.Trial_4.SLA),"--")
xline(height(SLA_timeseries) - height(SLA_Calculated.Retention.SLA),"--")
yline(0,"--")
ylim([-1,1])
title("Step Length Asymmetry Over Time With Visual Feedback")
ylabel("Step Length Asymmetry")
xlabel("Number of Strides")
hold off 

nexttile
plot(SB_timeseries);
hold on
xline(height(SLA_Calculated.Baseline.SB),"--")
xline(height(SLA_Calculated.Baseline.SB) + height(SLA_Calculated.Trial_1.SB),"--")
xline(height(SLA_Calculated.Baseline.SB) + height(SLA_Calculated.Trial_1.SB) + height(SLA_Calculated.Trial_2.SB),"--")
xline(height(SLA_Calculated.Baseline.SB) + height(SLA_Calculated.Trial_1.SB) + height(SLA_Calculated.Trial_2.SB) ...
    + height(SLA_Calculated.Trial_3.SB),"--")
xline(height(SLA_Calculated.Baseline.SB) + height(SLA_Calculated.Trial_1.SB) + height(SLA_Calculated.Trial_2.SB) ...
    + height(SLA_Calculated.Trial_3.SB) + height(SLA_Calculated.Trial_4.SB),"--")
xline(height(SB_timeseries) - height(SLA_Calculated.Retention.SB),"--")
yline(0,"--")
ylim([-1,1])
title("Step Length Asymmetry Over Time During Split Belt")
ylabel("Step Length Asymmetry")
xlabel("Number of Strides")
hold off 


% Plot Trajectory Data for SB
figure
tiledlayout (2,2)
sgtitle([SubjectID,' SB'])
for i =1:4
nexttile
plot(QTMFilt.SB{1,i}.Trajectories.X(:,7)) % R ankle
hold on
plot(QTMFilt.SB{1,i}.Trajectories.X(:,8)) % L ankle
% Add kinematic heel strikes
plot((GaitEvents{1,i}.SB.Kinematics.R_HS),QTMFilt.SB{1,i}.Trajectories.X(GaitEvents{1,i}.SB.Kinematics.R_HS,7),'x','Color','b')
%plot((GaitEvents{1,i}.SB.Kinetics.R_HS),QTMFilt.SB{1,i}.Trajectories.X(GaitEvents{1,i}.SB.Kinetics.R_HS,7),'x','Color','k')
plot((GaitEvents{1,i}.SB.Kinematics.L_HS),QTMFilt.SB{1,i}.Trajectories.X(GaitEvents{1,i}.SB.Kinematics.L_HS,8),'x', 'Color', [1, 0.5, 0])
%plot((GaitEvents{1,i}.SB.Kinetics.L_HS),QTMFilt.SB{1,i}.Trajectories.X(GaitEvents{1,i}.SB.Kinetics.L_HS,7),'x','Color','r')
xlim([1,2000])
title(['Trial ',num2str(i)])
ylabel('Ankle Trajectories')
xlabel('Sampling Over Time')
end 
legend ('R ankle', 'L ankle','R Kinematic Heel Strike','L Kinematic Heel Strike','Location','southeast outside')
hold off

% Plot Trajectory Data for SLA
figure
tiledlayout (2,2)
sgtitle([SubjectID,' SLA'])
for i =1:4
nexttile
plot(QTMFilt.SLA{1,i}.Trajectories.X(:,7)) % R ankle
hold on
plot(QTMFilt.SLA{1,i}.Trajectories.X(:,8)) % L ankle
% Add kinematic heel strikes
plot((GaitEvents{1,i}.SLA.Kinematics.R_HS),QTMFilt.SLA{1,i}.Trajectories.X(GaitEvents{1,i}.SLA.Kinematics.R_HS,7),'x','Color','b')
plot((GaitEvents{1,i}.SLA.Kinematics.L_HS),QTMFilt.SLA{1,i}.Trajectories.X(GaitEvents{1,i}.SLA.Kinematics.L_HS,8),'x', 'Color', [1, 0.5, 0])
xlim([1,2000])
title(['Trial ',num2str(i)])
ylabel('Ankle Trajectories')
xlabel('Sampling Over Time')
end 
legend ('R ankle', 'L ankle','R Kinematic Heel Strike','L Kinematic Heel Strike','Location','southeast outside')
hold off

% Plot Example Load Cell Data for SB
 % QTMFilt.SB{1,i).Load: the 9th column is RLoad and 10th column is LLoad
figure
tiledlayout (2,2)
sgtitle([SubjectID,' SB'])
for i =1:4
nexttile
plot(QTMFilt.SB{1,i}.Load(:,9)) % R ankle
hold on
plot(QTMFilt.SB{1,i}.Load(:,10)) % L ankle
% Add Kinetic Heel Strikes
plot((GaitEvents{1,i}.SB.Kinetics.R_HS),QTMFilt.SB{1,i}.Load(GaitEvents{1,i}.SB.Kinetics.R_HS,9),'x','Color','b')
plot((GaitEvents{1,i}.SB.Kinetics.L_HS),QTMFilt.SB{1,i}.Load(GaitEvents{1,i}.SB.Kinetics.L_HS,10),'x','Color','k')
title(['Trial ',num2str(i)])
ylabel('Load Cell Data')
xlabel('Sampling Over Time')
xlim([100000,200000])
end 
legend ('R ankle', 'L ankle', 'R Kinetic Heel Strike','L Kinetic Heel Strike','Location','southeast outside')
hold off

%  Plot Example Step Length Data for SB
figure
tiledlayout (2,2)
sgtitle([SubjectID,' SB'])
for i =1:4
nexttile
plot(StepLength{1,i}.SB.Right) % R ankle
hold on
plot(StepLength{1,i}.SB.Left) % L ankle
title(['Trial ',num2str(i)])
ylabel('Step Lengths')
xlabel('Strides')
end 
legend ('R ankle', 'L ankle')
hold off

%  Plot Example Step Length Data for SLA
figure
tiledlayout (2,2)
sgtitle([SubjectID,' SLA'])
for i =1:4
nexttile
plot(StepLength{1,i}.SLA.Right) % R ankle
hold on
plot(StepLength{1,i}.SLA.Left) % L ankle
title(['Trial ',num2str(i)])
ylabel('Step Lengths')
xlabel('Strides')
end 
legend ('R ankle', 'L ankle')
hold off