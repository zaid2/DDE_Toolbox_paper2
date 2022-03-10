function J = fcn_ctrl_dp(t,p)

% y = zeros(1,length(p),length(t));

% T0 = 1;
% T1 = t;
% T2 = 2*t.^2-1;
% 
% J = zeros(1,4,length(t));
% J(1,2,:) = T0;
% J(1,3,:) = T1;
% J(1,4,:) = T2;
% 
% % y = p(2).*T0+p(3).*T1+p(4).*T2;

J = zeros(1,size(p,1),length(t));
N = size(p,1);
T = zeros(N,length(t));
T(1,:) = ones(1,length(t));
T(2,:) = t;
for i=3:N
T(i,:) = 2*t.*T(i-1,:)-T(i-2,:);
end

for i=1:N

J(1,i,:) = T(i,:);

end

end







