function [data, y] = bc(prob, data, u)

T0 = u(1);
T  = u(2);
Omega = u(3); 

% 
y = [T0; T-2*pi/Omega];
% y = [T0; T-Tpar];


end