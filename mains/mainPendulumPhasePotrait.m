%%% Simple pendulum phase potrait
% Rishav (2020/9/7)
clc
clear
close all

% Simulation parameters
start_time = 0;
stop_time = 20;
dt = 0.1;
time  = start_time:dt:stop_time;

% Pendulum coefficients
model = 1; % 0 = linear and 1 = non linear
L = 5; % Length of pendulum
d = 0.5; % Damping coeff

for theta_init = -2*pi:0.5:2*pi
    for theta_dot_init = -3:0.5:3

% Initial state: [theta0;theta_dot0]
state = zeros(2,length(time));
state(:,1) = [theta_init;theta_dot_init];

% RK4 loop
for t = 1:length(time)-1
  if model == 1; fn = @(t,y)simplePendulum(t,y,L,d);
  elseif model == 0; fn = @(t,y)linearSimplePendulum(t,y,L,d); end
  state(:,t+1) = RK4(fn,state(:,t),dt,t);
end

plot(state(1,:),state(2,:)) % Phase plot
xlabel('theta'); ylabel('theta dot');
title('Simple Pendulum Phase Potrait');
hold on; drawnow;
    end
grid on;
end


