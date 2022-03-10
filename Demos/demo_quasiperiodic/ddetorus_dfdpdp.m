function J = ddetorus_dfdpdp(t, x, y, p)
%TORUS_DFDX   'coll'-compatible encoding of Jacobian of vector field w.r.t. problem variables

x1 = x(1,:);
x2 = x(2,:);

y1 = y(1,:);
y2 = y(2,:);

om1   = p(1,:);
Texc  = p(2,:);
alpha = p(3,:);

om2 = 2*pi./Texc;

r = sqrt(x1.^2+x2.^2);
J = zeros(2,3,3,numel(x1));

  J(1,2,2,:)= y1.*r.*(-cos(om2.*t)).*((-2*pi./Texc.^2).*t).*(-2*pi./Texc.^2).*t+...
              y1.*r.*(-sin(om2.*t)).*(-2*pi.*(-2./Texc.^3)).*t;
  J(2,2,2,:)= y2.*r.*(-cos(om2.*t)).*((-2*pi./Texc.^2).*t).*(-2*pi./Texc.^2).*t+...
              y2.*r.*(-sin(om2.*t)).*(-2*pi.*(-2./Texc.^3)).*t;
  
% y(1,:) = -om1.*x2+y1.*(1+r.*(cos(om2.*t)-1));
% y(2,:) = om1.*x1+y2.*(1+r.*(cos(om2.*t)-1));

end



