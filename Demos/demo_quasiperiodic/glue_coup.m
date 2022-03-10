function [data, y] = glue_coup(prob, data, u)


alpha     = u(end-1);
T = u(end);

y = zeros(data.nsegs,1);
for i=1:data.nsegs

eqn = u(i)-alpha/T;    

y(i,1) =  eqn;    
end


end


