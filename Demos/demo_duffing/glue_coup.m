function [data, y] = glue_coup(prob, data, u)

gamma1e = u(1);
T       = u(2);
alpha = 1;

y = gamma1e-alpha/T;
end


