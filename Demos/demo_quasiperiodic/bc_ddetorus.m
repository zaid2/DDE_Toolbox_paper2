function fbc = bc_ddetorus(data, T0, T, x0, x1, p)
% Copyright (C) Zaid Ahsan
%TORUS_BC   Torus boundary conditions.
%
% Trajectory end points lie on a curve on the invariant torus. The return
% map corresponds to identical times-of-flight and describes a rigid
% rotation.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: torus_bc.m 2839 2015-03-05 17:09:01Z fschild $

% yphase = int_phase(x0, data);
% fbc = [T0; T-2*pi/p(2); data.F*x1-data.RF*x0; yphase];
fbc = [T0; T-p(2); x0(2)];


end


function yphase = int_phase(x0,data)

x0ref = data.x0;
nsegs = data.nsegs;
dim = data.dim;

N = (nsegs-1)/2;   %====number of modes
NCOL=4;
[tc, wts] = coll_nodes(NCOL);

%====Calculation of x0 at the collocation nodes
Thj   = 2*pi*(0:1/(2*N+1):1); 
Th1 = Thj(1:end-1)+2*pi*(1+tc(:))/2/(2*N+1);
Th2 = kron(1:N, Th1(:));
F1inv = [ones((2*N+1)*NCOL,1) reshape([cos(Th2);sin(Th2)],...
                                                  [(2*N+1)*NCOL, 2*N])];
F2 = kron(F1inv,eye(dim)); 

coeff1 = data.F*x0;
x0_cn = F2*coeff1;

%====Calculation of x0ref at the collocation nodes
coeff2   = data.F*x0ref;
x0ref_cn = F2*coeff2;

%====Calculation of dx0ref at the collocation nodes
G1 = reshape([-sin(Th2);cos(Th2)],[(2*N+1)*NCOL, 2*N]);
G2 = repmat(kron((1:N), ones(1,2)), [(2*N+1)*NCOL,1]);
F3 = kron([zeros((2*N+1)*NCOL,1) G1.*G2], eye(dim));
x0refp_cn = F3*coeff2;

wts  = repmat(wts, [dim, 2*N+1]);
wts2 = spdiags(wts(:),0,(2*N+1)*NCOL*dim,(2*N+1)*NCOL*dim);
% Wts = kron(repmat(wts, [1,2*N+1]), ones(1,dim));

yphase = (x0_cn(:)-x0ref_cn(:))'*wts2*x0refp_cn(:)*2*pi/2/(2*N+1);

end


function [nds,wts]=coll_nodes(m)

n=(1:m-1)';
g=n.*(1./sqrt(4.*n.^2-1));
J=-diag(g,1)-diag(g,-1);

[w,x]=eig(J);
nds=diag(x);

wts=2*w(1,:).^2;
end






