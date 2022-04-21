function [t0,x0,y0] = init_guess_dde23(tau)

% Copyright (C) Zaid Ahsan

sol = dde23(@ddex1de,tau,@ddex1hist,[0, 2]);
% figure;
% plot(sol.x,sol.y(1,:))

t0 = (linspace(0,2))';
x0 = (deval(sol,t0))';

% Crude approximation of history
y0=zeros(length(t0),2);
for i=1:length(t0)
   if t0(i)<=1
      y0(i,1)=1;
   else
      y0(i,1)=interp1(t0,x0,t0(i)-1); 
   end
end

% Approximation of controls
p=[0,0,0];
y0(:,2) = p(1)+p(2).*t0+p(3).*(2*t0.^2-1);

end

function s = ddex1hist(t)
% Constant history function for DDEX1.
s = 1;
end

function dxdt = ddex1de(t,x,Z)
% Differential equations function for DDEX1.
xd = Z(1);

dxdt = t*x+xd;
end













