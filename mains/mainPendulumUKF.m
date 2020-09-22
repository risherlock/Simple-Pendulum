%%% Noisy simple pendulum tracking using EKF
% Rishav  (2020/9/2)
clc
clear
close all

%%% Simulation parameters
start_time = 0;
stop_time = 20;
dt = 0.04;
time  = start_time:dt:stop_time; 

%%% System parameters
g = 9.81; % Acceleration due to gravity
L = 1; % Length of pendulum
d = 0; % Damping coeff
init_state = [0,0.2]'; % [theta, theta_dot] 
x_prev = init_state; % Initial state
P = [0.001,0;0,0.001]; % Initial state covariance

%%% Define noise assumptions
Q = diag([0.01 0.01]); R = 0.01;
sys_params = {g,L,R,Q,d};
SNR = 25; 

% Initial state: [theta;theta_dot]
z = zeros(2,length(time));
state = zeros(2,length(time));
state(:,1) = init_state;
z(:,1) = [0.0;2.0];

%%% Generate ground truth
for t = 1:length(time)-1
  fn = @(t,y)simplePendulum(t,y,L,d);
  z(:,t+1) = RK4(fn,z(:,t),dt,t);
end
truth = z;
z = awgn(z,SNR); % Noisy measurement

%%% UKF
% UKF parameters initializations
n = 2; % Size of the state vector
alpha  = 1; % Primary scaling parameter
beta = 2; % Secondary scaline parameter
kappa = 0; % Tertiary scaling parameter

lambda = alpha^2*(n+kappa) - n;
W_mean = ones(2*n+1,1)*1/(2*(n+lambda));
W_cov = W_mean; W_mean(1) = lambda/(lambda+n);
W_cov(1) = lambda/(lambda+n) + 1 - alpha^2 + beta;
ukf_params = {W_mean,W_cov,lambda};

%%% Perform UKF
for i_iters  = 2:length(time)
    [state(:,i_iters),P] = pendulumUKF(state(:,i_iters-1),P,z(1,i_iters-1),dt,sys_params,ukf_params);
end

%%% Plot result
plot(time,truth(1,:),'LineWidth',1.5); hold on; % Ground truth
plot(time,z(1,:),'.','MarkerSize',10); hold on; % Measurement
plot(time,state(1,:),'LineWidth',1.5); grid on; % Estimation
legend('Ground truth','Measurement','Estimation');
title('Noisy Simple Pendulum UKF Estimation');