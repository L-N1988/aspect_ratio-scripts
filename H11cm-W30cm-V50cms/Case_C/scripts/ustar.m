clc; clear; close all;
casePath = '..';
load(fullfile(casePath, 'figure_data', 'u_stat.mat'));

%%
H = 11.02/100; % Water depth in m
W = 30/100; % Channel width in m
T = 14.9;  % Temperature in degree Celsius
S = 0.001435297; % Slope
log_index = Y > 0.02  & Y < 0.065;  % 0.0146 < y < 0.0517
% Linear regression to find the slope of the RSS profile
RSS_index = Y > 0.025 & Y < 0.08;

A = H * W; % Cross-sectional area
P = 2 * H + W; % Wetted perimeter
R = A / P; % Hydraulic radius

%% 
% The log layer velocity profile satisfies the following equation:
% U+ = 1/k * log(y+) + B
% where U+ = U/U_tau, y+ = y*U_tau/nu, k = 0.41, B = 5.2
% Use linear regression to find the best fit of U and y in the semilog plot.
U_log = U_xt(log_index);
Y_log = log(Y(log_index));
p = polyfit(Y_log, U_log, 1);
% The log layer velocity profile is plotted as a reference.
k = 0.41;
U_tau_log = k * p(1);
disp(['u_tau from log law fit: ', num2str(U_tau_log)]);
nu = 1e-6 * (60/(T+40))^1.45;  % Kinematic viscosity

%% 
figure();
scatter(U_xt, Y, 'o');
set(gca, 'YScale', 'log');  % Set log scale for y axis
xlabel("$U(\rm m/s)$", "FontSize", 14, 'Interpreter', 'latex');
ylabel("$y(\rm m)$", "FontSize", 14, 'Interpreter', 'latex');
set(gca, 'FontSize', 10); % Sets tick label size for both x and y axes
hold on;
% The log layer velocity profile is plotted as a reference.
plot(polyval(p, Y_log), Y(log_index), 'r--', 'LineWidth', 1.5);
legend('Data', 'Log law fit', 'Location', 'Best');

%%
% u_tau from slope of channel
% u_tau = sqrt(g * H * S)
u_tau_slope = sqrt(9.81 * R * S);
disp(['u_tau from slope of channel: ', num2str(u_tau_slope)]);

%%
% u_tau from RSS linear depth profile
figure();
scatter(RSS, Y, 'o');
xlabel("$-\overline{u'v'}(\rm m^2/s^2)$", "FontSize", 14, 'Interpreter', 'latex');
ylabel("$y(\rm m)$", "FontSize", 14, 'Interpreter', 'latex');
hold on;
RSS_linear = RSS(RSS_index);
Y_linear = Y(RSS_index);
p_trans = polyfit(Y_linear, RSS_linear, 1);
u_tau_RSS = sqrt(p_trans(2));
disp(['u_tau from RSS linear depth profile: ', num2str(u_tau_RSS)]);
% Add the linear fit to the plot
p = polyfit(RSS_linear, Y_linear, 1);
plot(RSS_linear, polyval(p, RSS_linear), 'r--', 'LineWidth', 1.5);
set(gca, 'FontSize', 10); % Sets tick label size for both x and y axes
legend('Data', 'RSS linear fit', 'Location', 'Best');