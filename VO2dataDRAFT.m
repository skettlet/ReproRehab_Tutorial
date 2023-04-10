% This code converts the measured VO2 and VCO2 data from an EXCEL file to a
% Matlab data file.

%% Load data
temp = dir('*.xls'); % Take .xls file and store it in a variable
%temp2 = xlsread(temp.name); % Read and store numeric values of the Excel file
temp2 = importdata(temp.name);
data = temp2.data;

%% Define constants
mass = data(6,9);
weight = mass*9.8;            % Weight in N
height = data(6,4)/100;       % Height in m
nrows = size(data,1);         % Extracting the number of rows in the matrix 'temp' which will be used to get the data values

%% Extract the actual data values from the Excel file
% To extract the data values from "temp", you should find the first and
% last row of the data. If you open "temp", you can see that the data are
% stored from the 29th to the 336th rows. 

% The following two lines of code find the end of the data by looking for
% NaN values.
na = find(isnan(data(29:nrows,1))); 
End_of_Data = na(1)+27;             % Last data point

% Store data values in a variable "VO2_data"
VO2_data = data(29:End_of_Data,:);


%% Create a vector containing all measured data including the trial type 
% Variable name in each column and corresponding units
% 1: Time(min) 
% 2: Heart_rate(not measured) 
% 3: VO2(L/min) 
% 4: VO2/kg(ml/(kg*min)) 
% 5: VCO2/kg(ml/(kg*min)) 
% 6: RER 
% 7: VE(L/min) 
% 8: Resting energy expenditure (kcal/day) 
% 9: FAT(%) 
% 10 CHO(%) 
% 11: Trial_type

time = VO2_data(:,1);
VO2Liters = VO2_data(:,2);
VO2Mill = VO2_data(:,3);
RER = VO2_data(:,7);

RER_smooth = smooth(RER);
VO2Mill_smooth = smooth(VO2Mill);
VO2Liters_smooth = smooth(VO2Liters);

figure
plot(time,RER)
hold on
plot(time,RER_smooth,"LineWidth",2)

figure
plot(time,VO2Liters)
hold on
plot(time,VO2Liters_smooth,"LineWidth",2)

figure
plot(time,VO2Mill)
hold on
plot(time,VO2Mill_smooth,"LineWidth",2)

%% Extract VO2 from walking and running trials

VO2_walking = VO2_data(VO2_data(:,11) == 2,[1 3]);
VO2_running = VO2_data(VO2_data(:,11) == 4,[1 3]);

Cals_Walking = Calories_Burned(VO2_walking(:,1),VO2_walking(:,2));
Cals_Running = Calories_Burned(VO2_running(:,1),VO2_running(:,2));

%--------------------------------------------------------------------------
function Total_Cals = Calories_Burned(Time,VO2)

Time_Intervals = diff(Time);
O2 = Time_Intervals.*(VO2(2:end));
Total_O2 = sum(O2);
Total_Cals = 4.9*Total_O2;
end