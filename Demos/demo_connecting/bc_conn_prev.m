function [data, y] = bc_conn_prev(prob, data, u)

x02 = u(1);
T0  = u(2);
x1  = u(3:4);
eta = u(5);

% wkplus = [-0.858781883599262; -0.512341367060774];
% lambdaplus = 0.596590795457267;

A0 = [0 1;1 eta];
[U,V]=eig(A0);
wkplus = U(:,2);
lambdaplus = V(2,2);

% wkplus = [-0.858642953046343; -0.512574169446585];
% lambdaplus = 0.596958453601750;

y = [x02; T0; -wkplus'*x1];

end