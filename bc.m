function [data, y] = bc(prob, data, u)

x0 = u(1);
T0 = u(2);

y = [x0-1; T0];

end