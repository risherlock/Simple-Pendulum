function [x_hat,P_hat] = pendulumUKF(x_prev,P_prev,z,dt,sys_params,ukf_params)
% Extended Kalman filter implementation for noisy pendulum tracking
%
% Input: n = 
% x = Present state i.e. previous estimate (nX1)
% u = Control input (px1)
% z = Measurement (mx1)
% sys_params = {F,B,H,P,Q,R}, a cell array 
%       F = State transition matrix (nXn)
%       B = Control matrix that maps control to state variables (nXp)
%       H = Measurement matrix that maps measurements to state (mXn)
%       P = State covariance matrix (n*n)
%       Q = Process covariance matrix (nX1)
%       R = Measurement covariance matrix (mXm)
% 
% Output:
%   state       = State estimate by KF (nX1)
%   covariance  = State covariance estimate by KF (nXn)
%
% Rishav (2020/9/3)
 
% Unpack system parameters
g = sys_params{1};
L = sys_params{2};
R = sys_params{3};
Q = sys_params{4};
d = sys_params{5};

% Unpack UKF params
W_mean = ukf_params{1};
W_cov = ukf_params{2};
lambda = ukf_params{3};
n = 2; m = 1;

%%% Prediction
sqrtP = chol(P_prev,'lower');
chi = [x_prev, ... % 2n+1 sigma points 
       x_prev + sqrt(n+lambda)*sqrtP, ...
       x_prev - sqrt(n+lambda)*sqrtP]; 

chi_x = zeros(n,2*n+1); 
chi_z = zeros(m,2*n+1);   
for i_iters = 1:2*n+1 
    % Extract states from sigma points
    theta = chi(1,i_iters);
    theta_dot = chi(2,i_iters);
    
    % Propagate sigma points through nonlinear process and observation model
    chi_x(:,i_iters) = pendulumPropagation([theta,theta_dot]',dt,g,L,d);
    chi_z(:,i_iters) = pendulumMeasurement(theta,L);
end
x_ = chi_x*W_mean; % Predicted state
z_ = chi_z*W_mean; % Predicted measurement

P_ = Q; Pyy = R; Pxy = zeros(n,m);  
for i_iters = 1:2*n+1 
    M = chi_x(:,i_iters) - x_;
    N = chi_z(:,i_iters) - z_;
   
   P_ = P_ + W_cov(i_iters)*(M*M'); % Predicted covariance
   Pyy = Pyy + W_cov(i_iters)*(N*N'); % Measurement/innovation covariance
   Pxy = Pxy + W_cov(i_iters)*(M*N'); % State-measurement cross-covariance
end

%%% Update
K = Pxy/Pyy; % Kalman gain
y_pre = z - z_; % Innovation or measurement pre-fit residual 
x_hat = x_ + K*y_pre; % Updated state estimate 
P_hat = P_ - K*Pyy*K'; % Updated estimate covariance 
end




