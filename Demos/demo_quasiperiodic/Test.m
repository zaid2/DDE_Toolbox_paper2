% Copyright (C) Zaid Ahsan
N=5;
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

check = real((kron(F,eye(2)))\(kron(R*F,eye(2))));



