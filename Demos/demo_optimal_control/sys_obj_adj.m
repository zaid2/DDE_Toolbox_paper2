function [data, J] = sys_obj_adj(prob, data,u)

NTST = data.ddaecoll.NTST;
NCOL = data.ddaecoll.NCOL;
W = data.W; 
xbp = u(data.xbp_idx);
ubp = u(data.xbp_idx(end)+1:end,1);

% ubp = ubp(data.tbp_idx);

xcn = W*xbp; 

xcndim = data.xcndim;
J = zeros(1,xcndim+NTST*(NCOL+1));

J(1,1:xcndim) = 2*xcn(:)';
J(1,xcndim+1:end) = 2*ubp(:)';
end