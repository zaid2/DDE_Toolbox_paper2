function Ac = Mat_coupling(p)
% Copyright (C) Zaid Ahsan
N=5; dim = 2;
varrho = 1/1.51111;
M = 2*N+1;
F = zeros(M,M);

jc = complex(0,1);
for i=1:M
    for j=1:M
      F(i,j) = exp(-2*pi*jc*(i-1)*(j-1)/M);  
    end
end

R = zeros(M,M);
R(1,1) = 1;
for i=1:N
R(i+1,i+1)     = exp(-2*pi*jc*i*varrho);
R(N+i+1,N+i+1) = exp(2*pi*jc*(N-(i-1))*varrho);
end

Ac = real((kron(F,eye(dim)))\(kron(R*F,eye(dim))));

% % varrho = p(3);
% varrho = 1/1.51111;
% N = 5;
% dim = 2;
% 
% Th = 2*pi*(0:2*N)/(2*N+1);
% Th = kron(1:N, Th');
% F  = [ones(2*N+1,1) 2*reshape([cos(Th);sin(Th)], [2*N+1 2*N])]'/(2*N+1);
% 
% Th  = (1:N)*2*pi*varrho;
% SIN = [ zeros(size(Th)) ; sin(Th) ];
% R   = diag([1 kron(cos(Th), [1 1])]);
% R   = R  + diag(SIN(:), +1)- diag(SIN(:), -1);
% 
% Th  = -(1:N)*2*pi*varrho;     %===difference between R and R1 is R1 has -2*pi*varrho
% SIN = [ zeros(size(Th)) ; sin(Th) ];
% R1   = diag([1 kron(cos(Th), [1 1])]);
% R1   = R1  + diag(SIN(:), +1)- diag(SIN(:), -1);
% 
% 
% 
% 
% Ac2 = inv(kron(F,eye(dim)))*kron(R1*F,eye(dim));
% check = Ac1-Ac2;
end
