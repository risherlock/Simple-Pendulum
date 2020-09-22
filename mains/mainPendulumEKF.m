%%% Noisy simple pendulum tracking using EKF
% Rishav  (2020/9/2)
clc
clear
close all

%%% Simulation parameters
start_time = 0;
stop_time = 10;
dt = 0.01;
time  = start_time:dt:stop_time; 

%%% System parameters
g = 9.81; % Acceleration due to gravity
L = 1; % Length of pendulum
d = 0.5; % Damping coeff
init_state = [0,0.2]'; % [theta, theta_dot] 
x_prev = init_state; % Initial state
P_prev = [0 ,0; 0, 0]; % Initial state covariance

%%% Define noise assumptions
Q = diag([0.001 0.001]); R = 2;
sys_params = {g,L,R,Q};
SNR = 25; 

% Initial state: [theta;theta_dot]
z = zeros(2,length(time));
state = zeros(2,length(time));
x_hat = zeros(2,length(time));
z(:,1) = [0.0;2.0];

%%% Generate ground truth
for t = 1:length(time)-1
  fn = @(t,y)simplePendulum(t,y,L,d);
  z(:,t+1) = RK4(fn,z(:,t),dt,t);
end
truth = z;
z = awgn(z,SNR); % Noisy measurement

%%% Perform EKF
for i_iters  = 2:length(time)
    [x_hat,P_hat] = pendulumEKF(x_prev,P_prev,z(1,i_iters),dt,sys_params);
    x_prev = x_hat; P_prev = P_hat;
    state(:,i_iters) = x_hat;
end

%%% Plot result
plot(time,truth(1,:),'LineWidth',1.5); hold on; % Ground truth
plot(time,z(1,:),'.','MarkerSize',10); hold on; % Measurement
plot(time,state(1,:),'LineWidth',1.5); grid on; % Estimation
legend('Ground truth','Measurement','Estimation');
title('Noisy Simple Pendulum EKF Estimation');