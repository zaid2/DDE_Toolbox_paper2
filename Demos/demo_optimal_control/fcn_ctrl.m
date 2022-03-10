function y = fcn_ctrl(t,p)


% T0 = 1;
% T1 = t;
% T2 = 2*t.^2-1;
% y = p(2).*T0+p(3).*T1+p(4).*T2;

N = size(p,1);
T = zeros(N,length(t));
T(1,:) = ones(1,length(t));
T(2,:) = t;
for i=3:N
T(i,:) = 2*t.*T(i-1,:)-T(i-2,:);
end
y = sum(p(1:end,:).*T,1);

end