function J = ddetorus_dfdxdp(t, x, y, p)
%TORUS_DFDX   'coll'-compatible encoding of Jacobian of vector field w.r.t. problem variables

x1 = x(1,:);
x2 = x(2,:);

y1 = y(1,:);
y2 = y(2,:);

om1  = p(1,:);
Texc = p(2,:);

om2  = 2*pi./Texc;

r = sqrt(x1.^2+x2.^2);

  J = zeros(2,2,3,numel(x1));

  J(1,2,1,:) = -1;
  J(2,1,1,:) =  1;
  
  J(1,1,2,:) = y1.*x1.*(-sin(om2.*t).*(-2*pi./Texc.^2).*t)./r;
  J(1,2,2,:) = y1.*(x2.*(-sin(om2.*t).*(-2*pi./Texc.^2).*t))./r;
  J(2,1,2,:) = y2.*x1.*(-sin(om2.*t).*(-2*pi./Texc.^2).*t)./r;
  J(2,2,2,:) = y2.*x2.*(-sin(om2.*t).*(-2*pi./Texc.^2).*t)./r;


end









