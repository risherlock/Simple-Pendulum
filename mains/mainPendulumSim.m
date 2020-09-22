%%% Simulation of the differential equation of the simple pendulum
% Rishav (2020/9/1)
clc
clear
close all

% Simulation parameters
start_time = 0;
stop_time = 25;
dt = 0.05;
time  = start_time:dt:stop_time; 

% Pendulum coefficients
model = 0; % 0 = linear and 1 = non linear
L = 3; % Length of pendulum
d = 0.5; % Damping coeff
theta_init = 0.1;
theta_dot_init = 0;
theta_l = 0;

% Initial state: [theta0;theta_dot0]
state = zeros(2,length(time));
state(:,1) = [theta_init;theta_dot_init];

% RK4 loop
for t = 1:length(time)-1
  if model == 1; fn = @(t,y)simplePendulum(t,y,L,d);
  elseif model == 0; fn = @(t,y)linearSimplePendulum(t,y,L,d); end
  state(:,t+1) = RK4(fn,state(:,t),dt,t);
end

% Plot
plot(time,state(1,:)); % Time plot
hold on; plot(time,state(2,:));
xlabel('Time'); 
legend('theta','theta dot');
title('Simple Pendulum');

figure;
plot(state(1,:),state(2,:)) % Phase plot
xlabel('Theta'); ylabel('Theta dot');
title('Simple Pendulum Phase Plot');

