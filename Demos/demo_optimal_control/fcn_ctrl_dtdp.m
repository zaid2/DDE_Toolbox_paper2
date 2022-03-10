function J = fcn_ctrl_dtdp(t,p)

% J = zeros(1,4,length(t));
% J(1,3,:) = 1;
% J(1,4,:) = 4*t;

% y = p(2).*T0+p(3).*T1+p(4).*T2;

J = zeros(1,size(p,1),length(t));

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

for i=1:N
  J(1,i,:) = Tp(i,:);  
end

end