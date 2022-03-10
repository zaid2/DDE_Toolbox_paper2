taucn_lshift = (Ti/Tj_k)*(taucn_k(:)-kron(ones(Ncn_k,1),Delta_k));
% mesh = 0:1/NTST:1;
% jcn_lshift = floor(interp1(mesh,1:NTST+1,taucn_lshift(:)));
rows = repmat((jcn_lshift(:)-1)'*(NCOL+1), [(NCOL+1),1]) +...
       repmat((1:(NCOL+1))',[1, length(taucn_lshift(:))]);
cols = repmat(1:length(taucn_lshift(:)),[NCOL+1,1]);
tc = 2*NTST*taucn_lshift(:,1)+1-2*jcn_lshift(:,1);
Lc = coll_L(seg.tm, tc); Lcp = coll_Lp(seg.tm, tc);
Lcn_lshift  = sparse(rows,cols,Lc',(NCOL+1)*NTST,length(taucn_lshift(:)));
Lcn_lshiftp = sparse(rows,cols,Lcp',(NCOL+1)*NTST,length(taucn_lshift(:)));

for s=1:Sk
    Aks  = Ak(:,1+(s-1)*dim:s*dim);
    Aksp = Akp(:,1+(s-1)*dim:s*dim,1:pdim);
    x_jks_bp = reshape(Xibp(1+(Xj_k(s)-1)*...
                        xbpdim:Xj_k(s)*xbpdim),[dim,NTST*(NCOL+1)]);
    
    J_T = J_T - (1/T)*Wmu_k'*Omega_k*kron(diag(taucn_lshift(:)...
                -Delta_k),eye(dim))*kron(Lcn_lshiftp',Aks)*x_jks_bp(:);
    
    J_Tjk = J_Tjk + (1/Tj_k)*Wmu_k'*Omega_k*...
                    kron(diag(taucn_lshift(:)-Delta_k),eye(dim))*...
                    kron(Lcn_lshiftp',Aks)*x_jks_bp(:);
    
    J_Deltak = J_Deltak+(Ti/Tj_k)*Wmu_k'*Omega_k*...
                            kron(Lcn_lshiftp',Aks)*x_jks_bp(:);
    
    J_p(:,1:pdim) = J_p(:,1:pdim)-Wmu_k'*Omega_k*...
                           kron(x_jks_bp*Lcn_lshift,eye(dim))'*...
                           reshape(Aksp(:),[dim^2,pdim]);
end