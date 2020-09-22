function [new_state] = pendulumPropagation(state,dt,g,L,d)
%%% Simple pendulum propagation function
% 2020/9/15

theta = state(1);
theta_dot = state(2);

new_state = [theta + theta_dot*dt, ...
            theta_dot - dt*g*sin(theta)/L - d*dt*theta_dot]';
end
 