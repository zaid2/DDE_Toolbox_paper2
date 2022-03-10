taucn_ushift = (Tj_k/Ti)*taucn_ks(:)+Delta_k;
mesh = 0:1/NTST:1;
jcn_ushift = floor(interp1(mesh,1:NTST+1,taucn_ushift(:)));

rows = repmat((jcn_ushift(:)-1)'*(NCOL+1), [(NCOL+1),1]) +...
    repmat((1:(NCOL+1))',[1, length(taucn_ushift(:))]);
cols = repmat(1:length(taucn_ushift(:)),[NCOL+1,1]);
tc = 2*NTST*taucn_ushift(:,1)+1-2*jcn_ushift(:,1);
Lc = coll_L(seg.tm, tc);
Lcn_ushift = sparse(rows,cols,Lc',(NCOL+1)*NTST,length(taucn_ushift(:)));

for s=1:Sk
    Aks = Ak(:,1+(s-1)*dim:s*dim);
    idx = 1+(idx_ks(1)-1)*dim+(Xj_k(s)-1)*xcndim:idx_ks(end)*dim+...
                                                     (Xj_k(s)-1)*xcndim;
    J_x(:,idx) = -(Tj_k/Ti)*kron(Lcn_ushift,Aks);
end