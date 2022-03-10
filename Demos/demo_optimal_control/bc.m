function [data, y] = bc(prob, data, u)

T0 = u(1);
T  = u(2); 

y = [T0; T-2];

end