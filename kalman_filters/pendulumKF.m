function [x_hat,P_hat] = pendulumKF(x_prev,P_prev,z,~,sys_params)
% Kalman filter implementation for noisy pendulum tracking
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
% Rishav (2020/9/14)
 
% Unpack system parameters
g = sys_params{1};
L = sys_params{2};
R = sys_params{3};
Q = sys_params{4};
d = sys_params{5};

% theta = x_prev(1);
theta = 0;
theta_dot = x_prev(2);
n = 2;

%% EKF implementation 
F = [0 1; -g*cos(theta)/L -d]; % State transition jacobian matrix 
H = [cos(theta), 0]; % Observation jacobian matrix 

% Prediction 
x_hat_  = F*[theta theta_dot]'; % (nX1)
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




