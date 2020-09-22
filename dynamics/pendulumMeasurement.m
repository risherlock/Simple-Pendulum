function [z] = pendulumMeasurement(theta,L)
%%% Simple pendulum propagation function
% 2020/9/15

z = L*sin(theta);
end