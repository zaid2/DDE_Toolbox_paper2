function init_guess_dde23()
% Copyright (C) Zaid Ahsan
alpha = 0.8255; T = 21.87; p1 = 0.5; p2 = -1.0783;
tau = alpha/T;
sol = dde23(@ddex1de,tau,@ddex1hist,[0, 1]);
figure(1)
plot(sol.y(1,:),sol.y(2,:))
hold on
plot(sol.y(1,end),sol.y(2,end),'ro')

% Stable manifold
phat2 = 0.5*(-p2+sqrt(4+p2^2));
v1 = phat2/(sqrt(phat2^2+1)); v2 = 1/(sqrt(phat2^2+1)); 
v = [v1;v2]; 
xs = -0.1:0.01:0.1;
ys = -(v2/v1)*xs;

hold on
plot(xs,ys)



end

function s = ddex1hist(t)
% Constant history function for DDEX1.

alpha = 0.8255; T = 21.87; p1 = 0.5; p2 = -1.0783;
eps = 5.51e-4;


phat1 = 0.5*(p2+sqrt(4+p2^2));
phat2 = 0.5*(-p2+sqrt(4+p2^2));
v1 = phat2/(sqrt(phat2^2+1)); v2 = 1/(sqrt(phat2^2+1)); l=phat1;

v = [v1;v2]; 

s = eps*v*exp(l*(T*t-alpha)); 

end

function dxdt = ddex1de(t,x,Z)

alpha = 0.8255; T = 21.87; p1 = 0.5; p2 = -1.0783;
dxdt = T*[x(2);x(1)-x(1)*Z(1,1)+p2*x(2)+p1*x(1)*x(2)];

end













