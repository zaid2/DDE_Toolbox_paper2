function J = fcn_ctrl_dt(t,p)

% Copyright (C) Zaid Ahsan

% T0 = 1;
% T1 = t;
% T2 = 2*t.^2-1;
% T3 = 4*t.^3-3*t;
% T4 = 8*t.^4-8*t.^2+1;
% T5 = 16*t.^5-20*t.^3+5*t;
% 
% J = zeros(1,length(t));
% J(1,:) = p(3,:)+p(4,:).*4.*t+p(5,:).*(12*t.^2-3)+p(6,:).*(32*t.^3-16*t)+p(7,:).*(80*t.^4);
% y = p(2).*T0+p(3).*T1+p(4).*T2;



N=size(p,1);

T = zeros(N,length(t));
T(1,:) = ones(1,length(t));
T(2,:) = t;
for i=3:N
T(i,:) = 2*t.*T(i-1,:)-T(i-2,:);
end

Tp = zeros(N,length(t));
Tp(1,:) = zeros(1,length(t));
Tp(2,:) = ones(1,length(t));
for i=3:N
Tp(i,:) = 2*T(i-1,:)+2*t.*Tp(i-1,:)-Tp(i-2,:);
end

J = sum(p(1:end,:).*Tp,1);
end





