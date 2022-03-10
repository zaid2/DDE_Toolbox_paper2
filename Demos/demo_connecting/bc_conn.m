function [data, y] = bc_conn(prob, data, u)

T0  = u(1);
x1  = u(2:3);
p2  = u(4);
lambda = u(5);
v1  =    u(6);
v2  =    u(7);
xbp_2dim =    u(8:end-1,1);
Tcr = u(end);

wk = [v1;v2];

%==phase condition
NTST = data.ddaecoll.NTST;
mesh = 0:1/NTST:1; 
m = 4;

L = zeros(1,NTST*(m+1));
k = floor(interp1(mesh,(1:NTST+1),Tcr));
% k = 69;
tc = 2*NTST*(Tcr-mesh(k))-1;

L(1,1+(m+1)*(k-1):(m+1)*k) = coll_L(data.tm,tc);



% wkplus = [-0.858642953046343; -0.512574169446585];
% lambdaplus = 0.596958453601750;

phat1 = 0.5*(p2+sqrt(4+p2^2));
phat2 = 0.5*(-p2+sqrt(4+p2^2));

% y = [T0; wk'*x1; lambda-phat1;...
%     v1-phat2/(sqrt(phat2^2+1)); v2-1/(sqrt(phat2^2+1)); L*xbp_2dim];

y = [T0; wk'*x1; lambda-phat1;...
    v1-phat2/(sqrt(phat2^2+1)); v2-1/(sqrt(phat2^2+1));L*xbp_2dim];

end