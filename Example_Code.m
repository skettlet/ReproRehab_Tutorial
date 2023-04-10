% This is a simple MATLAB code to calculate the average of two numbers and plot them

% Initialize the variables
num1 = 5;
num2 = 7;

% Calculate the average of the two numbers
sum = num1 + num2;
average = (sum) / 2;

% Display the result
disp(['The average of ' num2str(num1) ' and ' num2str(num2) ' is ' num2str(average)]);

% Plot the two numbers
figure;
bar([num1, num2]);
title('Numbers Plot');
xlabel('Numbers');
ylabel('Values');
