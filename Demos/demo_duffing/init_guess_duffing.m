function [t0,x0,y0] = init_guess_duffing(alpha,Omega,k,mu,zeta,a,b)

% Copyright (C) Zaid Ahsan

sol = dde23(@ddex1de,alpha,@ddex1hist,[0, 500],[],Omega,k,mu,zeta,a,b);
% figure;
% plot(sol.x,sol.y(1,:))

t   = (linspace(0,500,10000))';
sol0 = (deval(sol,t))';

T=2*pi/Omega;
ind_end = length(t);
val = t(end)-T;
[d,ind_begin] = min(abs(t-val));

x0 = sol0(ind_begin:ind_end,:);
t0 = t(ind_begin:ind_end,1);


% Crude approximation of history
y0=zeros(length(t0),2);
for i=1:length(t0)
    y0(i,1) = interp1(t,sol0(:,1),t0(i)-alpha);
    y0(i,2) = interp1(t,sol0(:,2),t0(i)-alpha);
end

t0 = t0-t0(1);
% x0 = x0';
% y0 = y0';
% figure;
% plot(t0/t0(end),x0(:,1),t0/t0(end),y0(:,1))
end

function s = ddex1hist(t,Omega,k,mu,zeta,a,b)
% Constant history function for DDEX1.
s = zeros(2,1);
end

function dxdt = ddex1de(t,x,Z,Omega,k,mu,zeta,a,b)
% Differential equations function for DDEX1.
x1lag = Z(1,1); x2lag = Z(2,1);
x1 = x(1); x2 = x(2);
dxdt = [x2;-2*zeta.*x2-x1-mu.*x1.^3+k.*x1lag+a.*cos(Omega.*t)+b.*sin(Omega.*t)];

end













