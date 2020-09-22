function state_dot = simplePendulum(~,state,L,d)
%%% Damped simple pendulum dynamics implementation
% Rishav (2020/9/2)

% Unpack state variables
theta = state(1);
theta_dot = state(2);
g = 9.8;

theta_dot_dot = -g*sin(theta)/L - d*theta_dot; 
state_dot = [theta_dot;theta_dot_dot];
end