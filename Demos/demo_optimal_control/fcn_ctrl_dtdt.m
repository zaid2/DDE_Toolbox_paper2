function J = fcn_ctrl_dtdt(t,p)

% Copyright (C) Zaid Ahsan

% J = zeros(1,length(t));
% J(1,:) = p(4,:)*4;

N = size(p,1);

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

Tpp = zeros(N,length(t));
Tpp(1,:) = zeros(1,length(t));
Tpp(2,:) = zeros(1,length(t));
for i=3:N
Tpp(i,:) = 2*Tp(i-1,:)+2*Tp(i-1,:)+2*t.*Tpp(i-1,:)-Tpp(i-2,:);    
end

J = sum(p(1:end,:).*Tpp,1);
end

























