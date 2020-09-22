function state_dot = linearSimplePendulum(~,state,L,d)
% Damped linear simple pendulum dynamics implementation
% The dynamics is linearized about theta_l
% Rishav (2020/9/7)

% Unpack state variables
theta = state(1);
g = 9.8;

%%% Non linear equations
% theta_dot_dot = -g*sin(theta)/L - d*theta_dot; 
% state_dot = [theta_dot;theta_dot_dot];

% Linearized state transition matrix
A = [0 1; -g*cos(theta)/L -d];
state_dot = A*state;
end