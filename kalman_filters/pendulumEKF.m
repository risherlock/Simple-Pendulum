function [x_hat,P_hat] = pendulumEKF(x_prev,P_prev,z,dt,sys_params)
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
% Rishav (2020/6/30)
 
% Unpack system parameters
g = sys_params{1};
L = sys_params{2};
R = sys_params{3};
Q = sys_params{4};
theta = x_prev(1);
theta_dot = x_prev(2);
n = 2;

%% EKF implementation 
F = [1 dt; -dt*g*cos(theta)/L 1]; % State transition jacobian matrix 
H = [cos(theta), 0]; % Observation jacobian matrix 

% Prediction 
x_hat_  = [theta+theta_dot*dt, theta_dot-dt*g*sin(theta)/L]'; % (nX1)
z_      = L*sin(theta); % Observation model           
P_      = F*P_prev*F' + Q; % Predicted covariance (nXn)

% Update 
S = H*P_*H' + R; % Innovation covariance (mX1)
C = P_*H'; % State-measurement cross-covariance (nxm)
K = C/S; % Optimal Kalman gain (nXm)

y_pre = z - z_; % Innovation or measurement pre-fit residual 
x_hat = x_hat_ + K*y_pre; % Updated state estimate 
P_hat = (eye(n) - K*H)*P_; % Updated estimate covariance 
% y_post  = z - H*x_hat; % Measurement post-fit residual 

end




