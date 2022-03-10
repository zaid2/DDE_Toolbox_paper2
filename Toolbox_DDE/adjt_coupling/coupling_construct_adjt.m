function [prob, data] = coupling_construct_adjt(prob, tbid, data, sol)
%COLL_CONSTRUCT_ADJT   Add COLL adjoint problem.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coll_construct_opt.m 2872 2015-08-06 20:12:06Z hdankowicz $

seg  = data.coupling_seg.ddaecoll_seg;
opt  = data.coupling_opt;

if ~isempty(sol.l0_coupling) 
   sol.l0_coupling = interp1(sol.tbp',sol.l0_coupling',seg.tbp','pchip')';
   sol.l0 = [sol.l0_coupling(:);sol.l0_bc(:)];
   if ~isempty(sol.tl0_coupling)
     sol.tl0_coupling = interp1(sol.tbp',sol.tl0_coupling',seg.tbp','pchip')';  
     sol.tl0 = [sol.tl0_coupling(:);sol.tl0_bc(:)];
   end
    
end

prob = coco_add_adjt(prob, tbid, @adj, @adj_DU, data, 'aidx', ...
    opt.guidx_adjt,...
    'l0', sol.l0(:), 'tl0', sol.tl0(:), 'adim', opt.adim);

if ~isempty(data.coupling_seg.pnames)
  pfid   = coco_get_id(tbid, 'pars');
  dnames = coco_get_id('d', data.coupling_seg.pnames);
  axidx  = coco_get_adjt_data(prob, tbid, 'axidx');
  prob   = coco_add_adjt(prob, pfid, dnames, 'aidx', axidx(opt.pnew_idx), ...
    'l0', sol.pars_l0, 'tl0', sol.pars_tl0);
end

end

%%
function [data, J] = adj(prob,data,u)

coupling_seg = data.coupling_seg;
coupling_opt = data.coupling_opt;
seg = coupling_seg.ddaecoll_seg;

NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
dim = seg.dim;
xbpdim = seg.xbpdim;
xcndim = seg.xcndim;

pdim = coupling_seg.pdim;

Xibp     = u(coupling_seg.Xibp_idx);
ybp      = u(coupling_seg.yibp_idx);
T        = u(coupling_seg.T_idx);
p        = u(coupling_seg.p_idx);
Delta    = u(coupling_seg.Delta_idx);
gamma    = u(coupling_seg.gamma_idx);

tau_pt = coupling_seg.tau_pt;

Ci = coupling_seg.Ci;

Xj_id     = coupling_seg.Xj_id;
Tj_id     = coupling_seg.Tj_id;
Ti_id     = coupling_seg.Ti_id;
Delta_id  = coupling_seg.Deltaj_id;
Nj        = coupling_seg.Nj;
xi        = coupling_seg.xi;

addim    = coupling_opt.addim;

taucn = seg.taucn;

%% Adjoint

if Nj~=0
    J_xcn   = sparse(xbpdim, Nj*xcndim);
else
    J_xcn = [];
end
J_T     = zeros(xbpdim,length(coupling_seg.T_idx));
J_p     = zeros(xbpdim,pdim);
J_Delta = zeros(xbpdim, length(coupling_seg.Delta_idx));
J_gamma = sparse(xbpdim, length(coupling_seg.gamma_idx));

kcn = floor(interp1([gamma(1:2:end-1);1],1:Ci+1,taucn));
kcn(kcn>Ci) = Ci;

for k=1:Ci
    
    Xj_k       = Xj_id{k};
    Sk         = length(Xj_k);
    Ti         = T(Ti_id);
    if Tj_id{k}~=0
    Tj_k       = T(Tj_id{k});
    end
    Delta_k    = Delta(Delta_id{k});
    
    gammab = gamma(2*k-1);
    gammae = gamma(2*k);
    
    idx_k      = find(kcn==k);
    taucn_k    = taucn(idx_k);
    Ncn_k      = length(taucn_k);
    
    idx      = 1+(idx_k(1)-1)*dim:idx_k(end)*dim;
    Omega_k  = seg.wts2(idx,idx);      %===Integration wts
    Wmu_k    = seg.W(idx,:);           %===Interpolation matrix for the Lagrange multiplier mu
    
    if Xj_k(1)==0
        t= Ti*taucn_k-Delta_k(1);
        dhistdt = coupling_seg.coupling{k}.func_dt(t(:)',repmat(p,[1,Ncn_k]));
        J_T(:,Ti_id) = J_T(:,Ti_id)-(1/2/NTST)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)),eye(dim))*dhistdt(:);
        J_Delta(:,Delta_id{k}(1)) = J_Delta(:,Delta_id{k}(1)) + (1/2/NTST)*Wmu_k'*Omega_k*dhistdt(:);
        
        dhistdp = coupling_seg.coupling{k}.func_dp(t(:)',repmat(p,[1,Ncn_k]));
        dprows = repmat(reshape(1:Ncn_k*dim, [dim Ncn_k]), [pdim 1]); % Index array for vectorization
        dpcols = repmat(1:pdim, [dim Ncn_k]);                         % Index array for vectorization
        dhistdp = sparse(dprows(:),dpcols(:),dhistdp(:));
        J_p(:,1:pdim) = J_p(:,1:pdim)-(1/2/NTST)*Wmu_k'*Omega_k*dhistdp;
    else
        
        if isempty(coupling_seg.coupling{k})
            Ak = eye(dim);
            Akp = zeros(dim,dim*Sk,pdim);
        else
            cfunc = coupling_seg.coupling{k}.func;
            cfunc_dp = coupling_seg.coupling{k}.func_dp;
            Mat    = cfunc(p);
            Mat_dp = cfunc_dp(p);
            
            rows  = coupling_seg.coupling{k}.rows;
            cols  = coupling_seg.coupling{k}.cols;
            
            if isempty(rows)
                Ak  = Mat;
                Akp = Mat_dp;
            else
                Ak  = Mat(rows,cols);
                Akp = Mat_dp(rows,cols,1:pdim);
            end
        end
        
        xi_bk = (Ti/Tj_k)*(gammab-Delta_k);
        xi_ek = (Ti/Tj_k)*(gammae-Delta_k);
            
        idx_ks   = intersect(find(taucn>=xi_bk),find(taucn<xi_ek));
        taucn_ks = taucn(idx_ks);
        
        if isempty(taucn_ks)
           continue;
        else
            taucn_ushift = (Tj_k/Ti)*taucn_ks(:)+Delta_k;
            mesh = 0:1/NTST:1;
            jcn_ushift = floor(interp1(mesh,1:NTST+1,taucn_ushift(:),'linear','extrap'));
            jcn_ushift(jcn_ushift>NTST) = NTST;
            jcn_ushift(jcn_ushift<=0)    = 1;
            
            rows = repmat((jcn_ushift(:)-1)'*(NCOL+1), [(NCOL+1),1]) +...
                repmat((1:(NCOL+1))',[1, length(taucn_ushift(:))]);
            cols = repmat(1:length(taucn_ushift(:)),[NCOL+1,1]);
            tc = 2*NTST*taucn_ushift(:,1)+1-2*jcn_ushift(:,1);
            Lc = coll_L(seg.tm, tc);
            Lcn_ushift = sparse(rows,cols,Lc',(NCOL+1)*NTST,length(taucn_ushift(:)));
                
            for s=1:Sk
                Aks  = Ak(:,1+(s-1)*dim:s*dim);
                idx = 1+(idx_ks(1)-1)*dim+(Xj_k(s)-1)*xcndim:idx_ks(end)*dim+(Xj_k(s)-1)*xcndim;
                J_xcn(:,idx) = -(Tj_k/Ti)*kron(Lcn_ushift,Aks);
            end
            
        end
        
        if isempty(taucn_k)
           continue;
        else
            %% delT, delTj, delDelta, delp
            taucn_lshift = (Ti/Tj_k)*(taucn_k(:)-kron(ones(Ncn_k,1),Delta_k));
            mesh = 0:1/NTST:1;
            jcn_lshift = floor(interp1(mesh,1:NTST+1,taucn_lshift,'linear','extrap'));   
            jcn_lshift(jcn_lshift>NTST) = NTST;
            jcn_lshift(jcn_lshift<=0)    = 1;
            
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
            x_jks_bp = reshape(Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim),[dim,NTST*(NCOL+1)]);
            
            J_T(:,Ti_id) = J_T(:,Ti_id) - (1/Ti)*Wmu_k'*...               %===Ti contribution
                Omega_k*kron(diag(taucn_lshift(:)),eye(dim))*...
                kron(Lcn_lshiftp',Aks)*x_jks_bp(:);
            
            J_T(:,Tj_id{k}) = J_T(:,Tj_id{k}) + (1/Tj_k)*Wmu_k'*...%===Tjks contribution
                Omega_k*kron(diag(taucn_lshift(:)),eye(dim))*...
                kron(Lcn_lshiftp',Aks)*x_jks_bp(:);
            
            J_Delta(:, Delta_id{k}) = J_Delta(:, Delta_id{k})+...
                                      (Ti/Tj_k)*Wmu_k'*Omega_k*...                    %===Deltaks contribution
                                       kron(Lcn_lshiftp',Aks)*x_jks_bp(:);
            
            J_p(:,1:pdim) = J_p(:,1:pdim)-Wmu_k'*Omega_k*kron(x_jks_bp*Lcn_lshift,eye(dim))'*reshape(Aksp(:),[dim^2,pdim]);
            end
                    
        end       
        
    end
end
     
     int_bcdim  = length(coupling_opt.X0_adjt_idx)+length(coupling_opt.X1_adjt_idx);
     J_coupling = [J_xcn eye(xbpdim) J_T sparse(xbpdim,int_bcdim) J_p J_Delta J_gamma];  %=== Adjoint of Coupling conditions + continuity of mu
     
%%  Adjoint contribution from Mesh conditions on gamma    
  J_fgamma = sparse(Ci+1,addim);
  gamma_adjt_idx = coupling_opt.gamma_adjt_idx;
for k=1:Ci+1
    if k==1
        J_fgamma(1,gamma_adjt_idx(1)) = 1;
    elseif 1<k && k<Ci+1
        %         gamma_eikminus1 = gamma(2*k-2); gamma_bik = gamma(2*k-1);
        %         fgamma(k) = gamma_eikminus1-gamma_bik;
        J_fgamma(k,gamma_adjt_idx(2*k-2)) =  1;
        J_fgamma(k,gamma_adjt_idx(2*k-1)) = -1;
    elseif k==Ci+1
        J_fgamma(k,gamma_adjt_idx(end)) = 1;
    end    
end


 %% Contributions from the internal boundary conditions
 J_int_bc = sparse((Ci-1)*dim,addim);
 
 for k=1:Ci-1
     Xj_k = Xj_id{k};
     
     if Xj_k(1)==0
         
         Xj_kplus1 = Xj_id{k+1};
         Skplus1 = length(Xj_kplus1);
         
         if isempty(coupling_seg.coupling{k+1})
             Akplus1  = eye(dim);
             Akplus1p = zeros(dim,Skplus1*dim,pdim);
         else
             cfunc = coupling_seg.coupling{k+1}.func;
             cfunc_dp = coupling_seg.coupling{k+1}.func_dp;
             Mat = cfunc(p);
             Mat_dp = cfunc_dp(p);
             rows =  coupling_seg.coupling{k+1}.rows;
             cols = coupling_seg.coupling{k+1}.cols;
             
             if isempty(rows)
                 Akplus1  = Mat;
                 Akplus1p = Mat_dp;
             else
                 Akplus1 = Mat(rows,cols); 
                 Akplus1p = Mat_dp(rows,cols,1:pdim);
             end
             
         end
         
         X0_adjt_idx = reshape(coupling_opt.X0_adjt_idx,[dim,Nj]);
         
         for s=1:Skplus1
             Akplus1s  = Akplus1(:,1+(s-1)*dim:s*dim);
             Akplus1sp = Akplus1p(:,1+(s-1)*dim:s*dim,:);
             
             % dfintbc_dxjkplus1sbp
             
             rows = repmat((1+(k-1)*dim:k*dim)', [1,dim]);
             cols = X0_adjt_idx(:,Xj_kplus1(s));
             cols = repmat(cols(:)',[dim,1]);
             J_int_bc = J_int_bc - sparse(rows(:),cols(:),Akplus1s, (Ci-1)*dim,addim);
             
             x_jkplus1s_bp = Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim);
             J_int_dp = -kron(coupling_seg.id0*x_jkplus1s_bp,eye(dim))'*reshape(Akplus1sp(:),[dim^2,pdim]);
             
             rows = repmat((1+(k-1)*dim:k*dim)', [1,dim*pdim]);
             cols = repmat(coupling_opt.p_adjt_idx(:)',[dim,1]);
             J_int_bc = J_int_bc  + sparse(rows,cols,J_int_dp,(Ci-1)*dim,addim);
         end
             
%              dhistdp = zeros(dim,pdim);
             dhistdp = coupling_seg.coupling{k}.func_dp(0,p);  
             rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim]);
             cols = repmat(coupling_opt.p_adjt_idx(:)',[dim,1]);
             J_int_bc = J_int_bc + sparse(rows,cols,dhistdp,(Ci-1)*dim,addim);
     else
         
         Xj_kplus1 = Xj_id{k+1};
         
         Sk      = length(Xj_k);
         Skplus1 = length(Xj_kplus1);
         
         X0_adjt_idx = reshape(coupling_opt.X0_adjt_idx,[dim,Nj]);
         X1_adjt_idx = reshape(coupling_opt.X1_adjt_idx,[dim,Nj]);
         
         if isempty(coupling_seg.coupling{k})
             Ak  = eye(dim);
             Akp = zeros(dim,dim*Sk,pdim);
         else
             cfunc = coupling_seg.coupling{k}.func;
             cfunc_dp = coupling_seg.coupling{k}.func_dp;
             Mat = cfunc(p);
             Mat_dp = cfunc_dp(p);
             
             rows =  coupling_seg.coupling{k}.rows;
             cols =  coupling_seg.coupling{k}.cols;
             
             if isempty(rows)
                 Ak  = Mat;
                 Akp = Mat_dp;
             else
                 Ak  = Mat(rows,cols);
                 Akp = Mat_dp(rows,cols,1:pdim);
%                  Akp = zeros(dim,dim*Sk,pdim);
             end  
         end
         
         if isempty(coupling_seg.coupling{k+1})
             Akplus1  = eye(dim);
             Akplus1p = zeros(dim,dim*Skplus1,pdim);
         else
             cfunc = coupling_seg.coupling{k+1}.func;
             cfunc_dp = coupling_seg.coupling{k+1}.func_dp;
             Mat = cfunc(p);
             Mat_dp = cfunc_dp(p);
             
             rows =  coupling_seg.coupling{k+1}.rows;
             cols =  coupling_seg.coupling{k+1}.cols;
             
             if isempty(rows)
                 Akplus1   = Mat;
                 Akplus1p = Mat_dp;
             else
                 Akplus1  = Mat(rows,cols);
                 Akplus1p = Mat_dp(rows,cols,1:pdim);
             end
         end
         
         % contribution w.r.t. x1_jks
         for s=1:Sk
             Aks  = Ak(:,1+(s-1)*dim:s*dim);
             Aksp = Akp(:,1+(s-1)*dim:s*dim,:);
             x_jks_bp = reshape(Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim),[dim,NTST*(NCOL+1)]);
             
             % contribution w.r.t. x1
             rows = repmat((1+(k-1)*dim:k*dim)', [1,dim]);
             cols = X1_adjt_idx(:,Xj_k(s));
             cols = repmat(cols(:)',[dim,1]);
             J_int_bc = J_int_bc + sparse(rows(:),cols(:),Aks, (Ci-1)*dim,addim);
             
             % contribution w.r.t. p
             rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim]);
             eqn_p = coupling_opt.p_adjt_idx;
             cols =  repmat(eqn_p(:)',[dim,1]);
             J_int_dp = kron(coupling_seg.id1*x_jks_bp(:),eye(dim))'*reshape(Aksp(:),[dim^2,pdim]);
             
             J_int_bc = J_int_bc + sparse(rows(:),cols(:),J_int_dp(:), (Ci-1)*dim,addim);
         end
         
         % contribution w.r.t. x0_jks
         for s=1:Skplus1
             Akplus1s  = Akplus1(:,1+(s-1)*dim:s*dim);
             Akplus1sp = Akplus1p(:,1+(s-1)*dim:s*dim,:);
             x_jkplus1s_bp = reshape(Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim),[dim,NTST*(NCOL+1)]);
             
             rows = repmat((1+(k-1)*dim:k*dim)', [1,dim]);
             cols = X0_adjt_idx(:,Xj_kplus1(s));
             cols = repmat(cols(:)',[dim,1]);
             J_int_bc = J_int_bc - sparse(rows(:),cols(:),Akplus1s, (Ci-1)*dim,addim);
             
             % contribution w.r.t. p
             rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim]);
             eqn_p = coupling_opt.p_adjt_idx;
             cols =  repmat(eqn_p(:)',[dim,1]);
             J_int_dp = -kron(coupling_seg.id0*x_jkplus1s_bp(:),eye(dim))'*reshape(Akplus1sp(:),[dim^2,pdim]);
             
             J_int_bc = J_int_bc + sparse(rows(:),cols(:),J_int_dp(:), (Ci-1)*dim,addim);
         end
     end
 end
 
 %% Adjoint contributions from mesh conditions on xi_bk and xi_ek
 J_xi_ek = [];
 J_xi_bk = [];
 for k=1:Ci-1
     Xj_k     = Xj_id{k};
    if Tj_id{k}~=0
    Tj_k     = T(Tj_id{k});
    end
    Delta_k  = Delta(Delta_id{k});
    
    Ti             = T(Ti_id);
    
    Xj_kplus1      = Xj_id{k+1}; 
    Tj_kplus1      = T(Tj_id{k+1});
    Delta_kplus1   = Delta(Delta_id{k+1});   
    
    if Xj_k(1)==0
        gamma_ek = gamma(2*k);
        dxieks_dTi = gamma_ek;
        dxieks_dgammaek = Ti;
        dxieks_dDeltaks  = -1;
        
        Ti_idx = coupling_opt.T_adjt_idx(Ti_id);
        gammaek_idx = coupling_opt.gamma_adjt_idx(2*k);
        Deltaks_idx  = coupling_opt.Delta_adjt_idx(Delta_id{k}(1));
        
        rows = ones(3,1); cols = [Ti_idx;gammaek_idx;Deltaks_idx];
        vals = [dxieks_dTi; dxieks_dgammaek; dxieks_dDeltaks];
        temp = sparse(rows,cols,vals,1,addim);
        J_xi_ek = [J_xi_ek;temp];
    else
        
        gamma_ek = gamma(2*k);
        dxiek_dTi = (1/Tj_k)*(gamma_ek-Delta_k);
        dxieik_dTjk = -(Ti/Tj_k^2)*(gamma_ek-Delta_k);
        dxiek_dgammaek = Ti/Tj_k;
        dxiek_dDeltak  = -(Ti/Tj_k);
        
        Ti_idx = coupling_opt.T_adjt_idx(Ti_id);
        Tjk_idx = coupling_opt.T_adjt_idx(Tj_id{k});
        gammaek_idx = coupling_opt.gamma_adjt_idx(2*k);
        Deltak_idx  = coupling_opt.Delta_adjt_idx(Delta_id{k});
        
        cols = [Ti_idx;Tjk_idx; gammaek_idx;Deltak_idx];
        rows = ones(4,1);
        vals = [dxiek_dTi dxieik_dTjk dxiek_dgammaek dxiek_dDeltak];
        temp = sparse(rows,cols,vals,1,addim);
        %             temp(1,idx) = [dxieks_dTi dxieiks_dTjks dxieks_dgammaek dxieks_dDeltaks];
        J_xi_ek = [J_xi_ek;temp];
        
    end
    
        gamma_bkplus1 = gamma(2*k+1);
        
        dxibkplus1_dTi            = (1/Tj_kplus1)*(gamma_bkplus1-Delta_kplus1);
        dxibkplus1_dTjkplus1      = -(Ti/Tj_kplus1^2)*(gamma_bkplus1-Delta_kplus1);
        dxibkplus1_dgammabkplus1  = (Ti/Tj_kplus1);
        dxibkplus1_dDeltakplus1   = -(Ti/Tj_kplus1);
    
        Ti_idx = coupling_opt.T_adjt_idx(Ti_id);
        Tjk_idx = coupling_opt.T_adjt_idx(Tj_id{k+1});
        gammabikplus1_idx = coupling_opt.gamma_adjt_idx(2*k+1);
        Deltakplus1_idx  = coupling_opt.Delta_adjt_idx(Delta_id{k+1});
        
        cols = [Ti_idx;Tjk_idx; gammabikplus1_idx;Deltakplus1_idx];
        rows = ones(4,1);
        vals = [dxibkplus1_dTi dxibkplus1_dTjkplus1 dxibkplus1_dgammabkplus1 dxibkplus1_dDeltakplus1];
        temp = sparse(rows,cols,vals,1,addim);
        %             temp(1,idx) = [dxibkplus1s_dTi dxibkplus1s_dTjkplus1s dxibkplus1s_dgammabkplus1 dxibkplus1s_dDeltakplus1s];
        J_xi_bk = [J_xi_bk;temp];
    
 end
 
 %% Overall adjoint problem
   J = [J_coupling;J_fgamma; J_int_bc; J_xi_ek; J_xi_bk];

end


%% Jacobian of the Adjoints
function [data,dJ] = adj_DU(prob,data,u)

coupling_seg = data.coupling_seg;
seg = coupling_seg.ddaecoll_seg;

opt = data.coupling_opt;

NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
dim = seg.dim;
xbpdim = seg.xbpdim;
xcndim = seg.xcndim;

pdim = coupling_seg.pdim;

Xibp    = u(coupling_seg.Xibp_idx);
ybp     = u(coupling_seg.yibp_idx);
T       = u(coupling_seg.T_idx);
p       = u(coupling_seg.p_idx);
Delta  = u(coupling_seg.Delta_idx);
gamma   = u(coupling_seg.gamma_idx);

tau_pt = coupling_seg.tau_pt;
taucn  = seg.taucn;

Ci = coupling_seg.Ci;

Xj_id     = coupling_seg.Xj_id;
Tj_id     = coupling_seg.Tj_id;
Ti_id     = coupling_seg.Ti_id;
Delta_id  = coupling_seg.Deltaj_id;
Nj        = coupling_seg.Nj;
xi        = coupling_seg.xi;

TDIM = length(coupling_seg.T_idx);
Deltadim = length(coupling_seg.Delta_idx);

addim = opt.addim;
dJ_coupling = sparse(xbpdim,addim*length(u));

kcn = floor(interp1([gamma(1:2:end-1);1],1:Ci+1,taucn));
kcn(kcn>Ci) = Ci;

for k=1:Ci
    Xj_k       = Xj_id{k};
    Sk         = length(Xj_k);
    Ti         = T(Ti_id);
    if Tj_id{k}~=0
        Tj_k       = T(Tj_id{k});
    end
    Delta_k    = Delta(Delta_id{k});
    
    gammab = gamma(2*k-1);
    gammae = gamma(2*k);
    
    idx_k      = find(kcn==k);
    taucn_k    = taucn(idx_k);
    Ncn_k      = length(taucn_k);
    
    idx      = 1+(idx_k(1)-1)*dim:idx_k(end)*dim;
    Omega_k  = seg.wts2(idx,idx);      %===Integration wts
    Wmu_k    = seg.W(idx,:);           %===Interpolation matrix for mu_cn
    
    if Xj_k(1)==0
        
        t  = Ti*taucn_k-Delta_k(1);
        dhistdtdt = coupling_seg.coupling{k}.func_dtdt(t(:)',repmat(p,[1,Ncn_k]));
        
        dhistdtdp = coupling_seg.coupling{k}.func_dtdp(t(:)',repmat(p,[1,Ncn_k]));

        dprows = repmat(reshape(1:Ncn_k*dim, [dim Ncn_k]), [pdim 1]); % Index array for vectorization
        dpcols = repmat(1:pdim, [dim Ncn_k]);                         % Index array for vectorization
        dhistdtdp = sparse(dprows(:),dpcols(:),dhistdtdp(:));
        
        dTi_dTi = -(1/2/NTST)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)),eye(dim))^2*dhistdtdt(:);
        dTi_dp  = -(1/2/NTST)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)),eye(dim))*dhistdtdp;
        dTi_dDelta  = +(1/2/NTST)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)),eye(dim))*dhistdtdt(:);
        
        dDelta_dDelta = -(1/2/NTST)*Wmu_k'*Omega_k*dhistdtdt(:);
        dDelta_dp = (1/2/NTST)*Wmu_k'*Omega_k*dhistdtdp;
        
        t  = Ti*taucn_k-Delta_k(1);
        dhistdpdp = coupling_seg.coupling{k}.func_dpdp(t(:)',repmat(p,[1,Ncn_k]));
  
        dpdprows = repmat(reshape(1:Ncn_k*dim, [dim Ncn_k]), [pdim^2 1]); % Index array for vectorization
        dpdpcols = repmat(1:pdim^2, [dim Ncn_k]);                         % Index array for vectorization
        dhistdpdp = sparse(dpdprows(:),dpdpcols(:),dhistdpdp(:));
        dp_dp = -(1/2/NTST)*Wmu_k'*Omega_k*dhistdpdp;
        
        % Jacobian of Ti contribution
        eqn_Ti = opt.T_adjt_idx(Ti_id);
        dTi_dTirows = 1:xbpdim;
        dTi_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Ti,[xbpdim,1]);
        dTi_dprows = repmat(1:xbpdim, [1,pdim]);
        dTi_dpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
            repmat(eqn_Ti,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
        dTi_dDeltarows = 1:xbpdim;
        dTi_dDeltacols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}(1)-1)*addim +...
            repmat(eqn_Ti,[xbpdim,1]);
        
        dTirows =  [dTi_dTirows(:);dTi_dDeltarows(:);dTi_dprows(:)];
        dTicols =  [dTi_dTicols(:);dTi_dDeltacols(:);dTi_dpcols(:)];
        dTivals =  [dTi_dTi(:);dTi_dDelta(:);dTi_dp(:)];
        
        dJ_coupling = dJ_coupling + sparse(dTirows(:),dTicols(:),dTivals,xbpdim,opt.dJcols);
        
        
        % Jacobian of Delta contribution
        eqn_Delta  = opt.Delta_adjt_idx(Delta_id{k}(1));
        dDelta_dTirows = 1:xbpdim;
        dDelta_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Delta,[xbpdim,1]);
        dDelta_dDeltarows = 1:xbpdim;
        dDelta_dDeltacols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}(1)-1)*addim + repmat(eqn_Delta,[xbpdim,1]);
        dDelta_dprows = repmat(1:xbpdim, [1,pdim]);
        dDelta_dpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
            repmat(eqn_Delta,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
        
        dDeltarows = [dDelta_dTirows(:);dDelta_dDeltarows(:);dDelta_dprows(:)];
        dDeltacols = [dDelta_dTicols(:);dDelta_dDeltacols(:);dDelta_dpcols(:)];
        dDeltavals = [dTi_dDelta(:);dDelta_dDelta(:);dDelta_dp(:)];
        
        dJ_coupling = dJ_coupling + sparse(dDeltarows(:),dDeltacols(:),dDeltavals,xbpdim,opt.dJcols);
        
        
        % Jacobian of p contribution
        eqn_p = opt.p_adjt_idx;
        dp_dTirows = repmat((1:xbpdim)',[1,pdim]);
        dp_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
        dp_dDeltarows = repmat((1:xbpdim)',[1,pdim]);
        dp_dDeltacols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}(1)-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
        dp_dprows = repmat((1:xbpdim)',[1,pdim*pdim]);
        dp_dpcols = (xbpdim*Nj+xbpdim+TDIM)*addim+...
            repmat(eqn_p(:)', xbpdim, pdim)+...
            repelem(0:addim:(pdim-1)*addim, xbpdim,pdim);
        dprows = [dp_dTirows(:);dp_dDeltarows(:);dp_dprows(:)];
        dpcols = [dp_dTicols(:);dp_dDeltacols(:);dp_dpcols(:)];
        dpvals = [dTi_dp(:);dDelta_dp(:);dp_dp(:)];
        
        dJ_coupling = dJ_coupling + sparse(dprows(:),dpcols(:),dpvals,xbpdim,opt.dJcols);
    else
        if isempty(coupling_seg.coupling{k})
            Ak = eye(dim);
            Akp = zeros(dim,dim*Sk,pdim);
            Akpp = zeros(dim,dim*Sk,pdim,pdim);
        else
            cfunc = coupling_seg.coupling{k}.func;
            cfunc_dp = coupling_seg.coupling{k}.func_dp;
            cfunc_dpdp = coupling_seg.coupling{k}.func_dpdp;
            
            Mat       = cfunc(p);
            Mat_dp    = cfunc_dp(p);
            Mat_dpdp  = cfunc_dpdp(p);
            rows =  coupling_seg.coupling{k}.rows;
            cols =  coupling_seg.coupling{k}.cols;
            
            if isempty(rows)
                Ak   = Mat;
                Akp  = Mat_dp;
                Akpp = Mat_dpdp;
            else
                Ak   = Mat(rows,cols);
                Akp  = Mat_dp(rows,cols,1:pdim);
                Akpp = Mat_dpdp(rows,cols,1:pdim,1:pdim);
%                 Akpp = zeros(dim,dim*Sk,pdim,pdim);
            end
            
%             Ak   = Mat(rows,cols);
%             Akp  = zeros(dim,dim*Sk,pdim);
%             Akpp = zeros(dim,dim*Sk,pdim,pdim);
        end
        
        
        %% Jacobian of Adjoint contribution associated with xj_ks
        xi_bk = (Ti/Tj_k)*(gammab-Delta_k);
        xi_ek = (Ti/Tj_k)*(gammae-Delta_k);
        
        idx_ks   = intersect(find(taucn>=xi_bk),find(taucn<xi_ek));
        taucn_ks = taucn(idx_ks);
        Ncn_ks   = length(taucn_ks);
        
        if isempty(taucn_ks)
            continue;
        else
            taucn_ushift = (Tj_k/Ti)*taucn_ks(:)+Delta_k;
%             jcn_ushift = shift_loc_nodes(taucn_ushift,tau_pt);
%             taucn_ushift(taucn_ushift<0) = 0; 
%             taucn_ushift(taucn_ushift>1) = 1;
            mesh = 0:1/NTST:1;
            jcn_ushift = floor(interp1(mesh,1:NTST+1,taucn_ushift(:),'linear','extrap'));
            jcn_ushift(jcn_ushift>NTST) = NTST;
            jcn_ushift(jcn_ushift<=0)   = 1;
            
            rows = repmat((jcn_ushift(:)-1)'*(NCOL+1), [(NCOL+1),1]) +...
                repmat((1:(NCOL+1))',[1, length(taucn_ushift(:))]);
            cols = repmat(1:length(taucn_ushift(:)),[NCOL+1,1]);
            tc = 2*NTST*taucn_ushift(:,1)+1-2*jcn_ushift(:,1);
            Lc =  coll_L(seg.tm, tc);
            Lcp = coll_Lp(seg.tm, tc);
            Lcn_ushift  = sparse(rows,cols,Lc',(NCOL+1)*NTST,length(taucn_ushift(:)));
            Lcn_ushiftp = sparse(rows,cols,Lcp',(NCOL+1)*NTST,length(taucn_ushift(:)));
            
            for s=1:Sk
                
                Aks   = Ak(:,1+(s-1)*dim:s*dim);
                Aksp  = Akp(:,1+(s-1)*dim:s*dim,1:pdim);
                
                dxjks_dTi      = (Tj_k/Ti^2)*kron(Lcn_ushift,Aks)+(Tj_k^2/Ti^3)*2*NTST*kron(Lcn_ushiftp*diag(taucn_ks(:)),Aks);
                dxjks_dTjk    = -(1/Ti)*kron(Lcn_ushift,Aks)-(Tj_k/Ti^2)*2*NTST*kron(Lcn_ushiftp*diag(taucn_ks(:)),Aks);
                dxjks_dDeltak = -(Tj_k/Ti)*2*NTST*kron(Lcn_ushiftp,Aks);
                
                dxjks_dp = [];
                for j=1:pdim
                    temp =  -(Tj_k/Ti)*kron(Lcn_ushift,Aksp(:,:,j));
                    dxjks_dp = [dxjks_dp;temp(:)];
                end
                
                rows = repmat((1:NTST*(NCOL+1)*dim)',[1,Ncn_ks*dim]);
                cols = repmat(1+(idx_ks(1)-1)*dim+(Xj_k(s)-1)*xcndim:idx_ks(end)*dim+(Xj_k(s)-1)*xcndim,[NTST*(NCOL+1)*dim,1]);
                
                dxjks_dTirows    = rows(:);
                dxjks_dTicols    = cols(:) + (xbpdim*Nj+xbpdim+Ti_id-1)*addim;
                dxjks_dTjkrows  = rows(:);
                dxjks_dTjkcols  = cols(:) + (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim;
                dxjks_dDeltakrows = rows(:);
                dxjks_dDeltakcols = cols(:) + (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}-1)*addim;
                dxjks_dprows = repmat(rows(:),[1,pdim]);
                dxjks_dpcols = repmat(cols(:),[1,pdim])+(xbpdim*Nj+xbpdim+TDIM)*addim+repmat((0:1:pdim-1)*addim,[xbpdim*Ncn_ks*dim,1]);
                
                dxjksrows = [dxjks_dTirows(:);dxjks_dTjkrows(:); dxjks_dDeltakrows(:); dxjks_dprows(:)];
                dxjkscols = [dxjks_dTicols(:);dxjks_dTjkcols(:); dxjks_dDeltakcols(:); dxjks_dpcols(:)];
                dxjksvals = [dxjks_dTi(:);dxjks_dTjk(:); dxjks_dDeltak(:); dxjks_dp(:)];
                
                dJ_coupling = dJ_coupling + sparse(dxjksrows,dxjkscols,dxjksvals,xbpdim,opt.dJcols);
            end
            
            
            if isempty(taucn_k)
                continue;
            else
                
                taucn_lshift = (Ti/Tj_k)*(taucn_k(:)-kron(ones(Ncn_k,1),Delta_k));
                jcn_lshift   = shift_loc_nodes(taucn_lshift,tau_pt);
                
                rows = repmat((jcn_lshift(:)-1)'*(NCOL+1), [(NCOL+1),1]) +...
                    repmat((1:(NCOL+1))',[1, length(taucn_lshift(:))]);
                cols = repmat(1:length(taucn_lshift(:)),[NCOL+1,1]);
                tc = 2*NTST*taucn_lshift(:,1)+1-2*jcn_lshift(:,1);
                Lc = coll_L(seg.tm, tc);
                Lcp = coll_Lp(seg.tm, tc);
                Lcpp = coll_Lpp(seg.tm, tc);
                Lcn_lshift   = sparse(rows,cols,Lc',(NCOL+1)*NTST,length(taucn_lshift(:)));
                Lcn_lshiftp  = sparse(rows,cols,Lcp',(NCOL+1)*NTST,length(taucn_lshift(:)));
                Lcn_lshiftpp = sparse(rows,cols,Lcpp',(NCOL+1)*NTST,length(taucn_lshift(:)));
                
                for s=1:Sk
                    Aks   = Ak(:,1+(s-1)*dim:s*dim);
                    Aksp  = Akp(:,1+(s-1)*dim:s*dim,1:pdim);
                    Akspp = Akpp(:,1+(s-1)*dim:s*dim,1:pdim,1:pdim);
                    
                    x_jks_bp = reshape(Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim),[dim,NTST*(NCOL+1)]);
                    
                    %% ===Jacobian of Tj_ks contribution
                    eqn_Tjk = opt.T_adjt_idx(Tj_id{k});
                    
                    dTjk_dxjksbp = (Ti/Tj_k^2)*Wmu_k'*...
                        Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        kron(Lcn_lshiftp',Aks);
                    
                    dTjk_dTi = (1/Tj_k^2)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        kron(Lcn_lshiftp',Aks)*x_jks_bp(:) +...
                        (Ti/Tj_k^3)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        2*NTST*kron((Lcn_lshiftpp*diag(taucn_k(:)-Delta_k))',Aks)*x_jks_bp(:);
                    
                    dTjk_dTjk = (-2*Ti/Tj_k^3)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        kron(Lcn_lshiftp',Aks)*x_jks_bp(:)...
                        -(Ti^2/Tj_k^4)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        2*NTST*kron((Lcn_lshiftpp*diag(taucn_k(:)-Delta_k))',Aks)*x_jks_bp(:);
                    
                    dTjk_dDeltak = -(Ti^2/Tj_k^3)*Wmu_k'*...
                        Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        2*NTST*kron(Lcn_lshiftpp',Aks)*x_jks_bp(:)-(Ti/Tj_k^2)*Wmu_k'*...%===Tjks contribution
                                                           Omega_k*kron(Lcn_lshiftp',Aks)*x_jks_bp(:);
                    
                    dTjk_dp = (Ti/Tj_k^2)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        kron(x_jks_bp*Lcn_lshiftp,eye(dim))'*reshape(Aksp(:),[dim^2,pdim]);
                    
                    dTjk_dxjksbprows = repmat((1:xbpdim)',[1,xbpdim]);
                    dTjk_dxjksbpcols = repmat(eqn_Tjk,[xbpdim,xbpdim])+...
                        (Xj_k(s)-1)*xbpdim*addim+repmat(0:addim:(xbpdim-1)*addim,[xbpdim,1]);
                    
                    dTjk_dTirows = 1:xbpdim;
                    dTjk_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Tjk,[xbpdim,1]);
                    dTjk_dTjkrows = 1:xbpdim;
                    dTjk_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + repmat(eqn_Tjk,[xbpdim,1]);
                    dTjk_dDeltakrows = 1:xbpdim;
                    dTjk_dDeltakcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}-1)*addim + repmat(eqn_Tjk,[xbpdim,1]);
                    dTjk_dprows = repmat(1:xbpdim, [1,pdim]);
                    dTjk_dpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
                        repmat(eqn_Tjk,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
                    
                    dTjkvals = [dTjk_dxjksbp(:);dTjk_dTi(:);dTjk_dTjk(:);dTjk_dDeltak(:);dTjk_dp(:)];
                    dTjkrows = [dTjk_dxjksbprows(:);dTjk_dTirows(:);dTjk_dTjkrows(:);dTjk_dDeltakrows(:);dTjk_dprows(:)];
                    dTjkcols = [dTjk_dxjksbpcols(:);dTjk_dTicols(:);dTjk_dTjkcols(:);dTjk_dDeltakcols(:);dTjk_dpcols(:)];
                    
                    dJ_coupling = dJ_coupling + sparse(dTjkrows,dTjkcols,dTjkvals,xbpdim,opt.dJcols);
                    
                    %% Jacobian of Delta_ks contribution
                    eqn_Deltak = opt.Delta_adjt_idx(Delta_id{k});
                    
                    dDeltak_dxjksbp = (Ti/Tj_k)*Wmu_k'*Omega_k* kron(Lcn_lshiftp',Aks);
                    dDeltak_dTi = (1/Tj_k)*Wmu_k'*Omega_k*kron(Lcn_lshiftp',Aks)*x_jks_bp(:) + (Ti/Tj_k^2)*Wmu_k'*Omega_k*...
                        2*NTST*kron((Lcn_lshiftpp*diag(taucn_k(:)-Delta_k))',Aks)*x_jks_bp(:);
                    dDeltak_dTjk = dTjk_dDeltak;
                    dDeltak_dDeltak = -(Ti^2/Tj_k^2)*Wmu_k'*Omega_k*2*NTST*kron(Lcn_lshiftpp',Aks)*x_jks_bp(:);
                    dDeltak_dp = (Ti/Tj_k)*Wmu_k'*Omega_k*kron(x_jks_bp*Lcn_lshiftp,eye(dim))'*reshape(Aksp(:),[dim^2,pdim]);
                    
                    dDeltak_dxjksbprows = repmat((1:xbpdim)',[1,xbpdim]);
                    dDeltak_dxjksbpcols = repmat(eqn_Deltak,[xbpdim,xbpdim])+...
                        (Xj_k(s)-1)*xbpdim*addim+repmat(0:addim:(xbpdim-1)*addim,[xbpdim,1]);
                    dDeltak_dTirows = 1:xbpdim;
                    dDeltak_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Deltak,[xbpdim,1]);
                    dDeltak_dTjkrows = 1:xbpdim;
                    dDeltak_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + repmat(eqn_Deltak,[xbpdim,1]);
                    dDeltak_dDeltakrows = 1:xbpdim;
                    dDeltak_dDeltakcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}-1)*addim + repmat(eqn_Deltak,[xbpdim,1]);
                    dDeltak_dprows = repmat(1:xbpdim, [1,pdim]);
                    dDeltak_dpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
                        repmat(eqn_Deltak,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
                    
                    dDeltakvals = [dDeltak_dxjksbp(:);dDeltak_dTi(:);dDeltak_dTjk(:);dDeltak_dDeltak(:);dDeltak_dp(:)];
                    dDeltakrows = [dDeltak_dxjksbprows(:);dDeltak_dTirows(:);dDeltak_dTjkrows(:);dDeltak_dDeltakrows(:);dDeltak_dprows(:)];
                    dDeltakcols = [dDeltak_dxjksbpcols(:);dDeltak_dTicols(:);dDeltak_dTjkcols(:);dDeltak_dDeltakcols(:);dDeltak_dpcols(:)];
                    
                    dJ_coupling = dJ_coupling + sparse(dDeltakrows,dDeltakcols,dDeltakvals,xbpdim,opt.dJcols);
                    
                    %% Jacobian of Ti contribution
                    eqn_Ti = opt.T_adjt_idx(Ti_id);
                    
                    dTi_dxjksbp  = -(1/Tj_k)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        kron(Lcn_lshiftp',Aks);
                    dTi_dTjk    = dTjk_dTi;
                    dTi_dDeltak = dDeltak_dTi;
                    dTi_dTi      = -(1/Tj_k^2)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        2*NTST*kron((Lcn_lshiftpp*diag(taucn_k(:)-Delta_k))',Aks)*x_jks_bp(:);
                    dTi_dp =       -(1/Tj_k)*Wmu_k'*Omega_k*kron(diag(taucn_k(:)-Delta_k),eye(dim))*...
                        kron(x_jks_bp*Lcn_lshiftp,eye(dim))'*reshape(Aksp(:),[dim^2,pdim]);
                    
                    dTi_dxjksbprows = repmat((1:xbpdim)',[1,xbpdim]);
                    dTi_dxjksbpcols = repmat(eqn_Ti,[xbpdim,xbpdim])+...
                        (Xj_k(s)-1)*xbpdim*addim+repmat(0:addim:(xbpdim-1)*addim,[xbpdim,1]);
                    dTi_dTirows = 1:xbpdim;
                    dTi_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Ti,[xbpdim,1]);
                    dTi_dTjkrows = 1:xbpdim;
                    dTi_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + repmat(eqn_Ti,[xbpdim,1]);
                    dTi_dDeltakrows = 1:xbpdim;
                    dTi_dDeltakcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}-1)*addim +...
                        repmat(eqn_Ti,[xbpdim,1]);
                    dTi_dprows = repmat(1:xbpdim, [1,pdim]);
                    dTi_dpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
                        repmat(eqn_Ti,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
                    
                    dTivals = [dTi_dxjksbp(:);dTi_dTi(:);dTi_dTjk(:);dTi_dDeltak(:);dTi_dp(:)];
                    dTirows = [dTi_dxjksbprows(:);dTi_dTirows(:);dTi_dTjkrows(:);dTi_dDeltakrows(:);dTi_dprows(:)];
                    dTicols = [dTi_dxjksbpcols(:);dTi_dTicols(:);dTi_dTjkcols(:);dTi_dDeltakcols(:);dTi_dpcols(:)];
                    
                    dJ_coupling = dJ_coupling + sparse(dTirows(:),dTicols(:),dTivals,xbpdim,opt.dJcols);
                    
                    %% Jacobian of p contribution
                    eqn_p = opt.p_adjt_idx;
                    
                    dp_dxjksbp = [];
                    for i=1:pdim
                        temp = -Wmu_k'*Omega_k*kron(Lcn_lshift',Aksp(:,:,i));
                        dp_dxjksbp = [dp_dxjksbp;temp(:)];
                    end
                    
                    dp_dTi = dTi_dp;
                    dp_dTjk = dTjk_dp;
                    dp_dDeltak = dDeltak_dp;
                    dp_dp = -Wmu_k'*Omega_k*kron(x_jks_bp*Lcn_lshift,eye(dim))'*reshape(Akspp(:),[dim^2,pdim*pdim]);
                    
                    rows = repmat((1:xbpdim)',[1,xbpdim]);
                    dp_dxjksbprows = repmat(rows,[1,pdim]);
                    cols = repmat(1:xbpdim,[xbpdim,1])+repmat(0:addim:(xbpdim-1)*addim,[xbpdim,1]);
                    dp_dxjksbpcols = repmat(cols(:),[1,pdim]) + repmat(eqn_p(:)'-1,[xbpdim*xbpdim,1]);
                    
                    dp_dTirows = repmat((1:xbpdim)',[1,pdim]);
                    dp_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
                    dp_dTjkrows = repmat((1:xbpdim)',[1,pdim]);
                    dp_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
                    dp_dDeltakrows = repmat((1:xbpdim)',[1,pdim]);
                    dp_dDeltakcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
                    dp_dprows = repmat((1:xbpdim)',[1,pdim*pdim]);
                    dp_dpcols = (xbpdim*Nj+xbpdim+TDIM)*addim+...
                        repmat(eqn_p(:)', xbpdim, pdim)+...
                        repelem(0:addim:(pdim-1)*addim, xbpdim,pdim);
                    
                    dpvals = [dp_dxjksbp(:);dp_dTi(:);dp_dTjk(:);dp_dDeltak(:);dp_dp(:)];
                    dprows = [dp_dxjksbprows(:);dp_dTirows(:);dp_dTjkrows(:);dp_dDeltakrows(:);dp_dprows(:)];
                    dpcols = [dp_dxjksbpcols(:);dp_dTicols(:);dp_dTjkcols(:);dp_dDeltakcols(:);dp_dpcols(:)];
                    
                    dJ_coupling = dJ_coupling + sparse(dprows(:),dpcols(:),dpvals,xbpdim,opt.dJcols);
                    
                end
                
            end
        end
        
        
        
    end
end


%% Jacobian of Adjoint conntribution from Mesh conditions on gamma
dJ_fgamma = sparse(Ci+1,opt.dJcols);

%% Jacobian of contributions from internal boundary conditions
dJ_int_bc = sparse((Ci-1)*dim,addim*length(u));
for k=1:Ci-1
    Xj_k = Xj_id{k};
    
    if Xj_k(1)==0
         Xj_kplus1 = Xj_id{k+1};
         Skplus1 = length(Xj_kplus1);
         
        if isempty(coupling_seg.coupling{k+1})
            Akplus1  = eye(dim);
            Akplus1p = zeros(dim,dim*Skplus1,pdim);
            Akplus1pp = zeros(dim,dim*Skplus1,pdim,pdim);
        else
            cfunc = coupling_seg.coupling{k+1}.func;
            cfunc_dp   = coupling_seg.coupling{k+1}.func_dp;
            cfunc_dpdp = coupling_seg.coupling{k+1}.func_dpdp;
            
            Mat      = cfunc(p);
            Mat_dp   = cfunc_dp(p);
            Mat_dpdp = cfunc_dpdp(p);
            
            rows =  coupling_seg.coupling{k+1}.rows;
            cols =  coupling_seg.coupling{k+1}.cols;
            
            if isempty(rows)
            Akplus1   = Mat;
            Akplus1p  = Mat_dp;
            Akplus1pp = Mat_dpdp;
            else
            Akplus1     = Mat(rows,cols);
            Akplus1p    = Mat_dp(rows,cols,1:pdim);
            Akplus1pp   = Mat_dpdp(rows,cols,1:pdim,1:pdim);
            end
            

        end
         
        X0_adjt_idx = reshape(opt.X0_adjt_idx,[dim,Nj]); 
        
                % contribution w.r.t. x0_jks
        for s=1:Skplus1
            Akplus1s   = Akplus1(:,1+(s-1)*dim:s*dim);
            Akplus1sp  = Akplus1p(:,1+(s-1)*dim:s*dim,:);
            Akplus1spp = Akplus1pp(:,1+(s-1)*dim:s*dim,:,:);
            x_jkplus1s_bp = reshape(Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim),[dim,NTST*(NCOL+1)]);
            
            % Jacobian w.r.t. p
            rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim*pdim]);
            eqn_p = opt.p_adjt_idx;
            cols = (xbpdim*Nj+xbpdim+TDIM)*addim+repmat(eqn_p(:)', dim, pdim)+...
                repelem(0:addim:(pdim-1)*addim, dim,pdim);
            dJ_int_dp = -kron(coupling_seg.id0*x_jkplus1s_bp(:),eye(dim))'*reshape(Akplus1spp(:),[dim^2,pdim*pdim]);
            
            dJ_int_bc = dJ_int_bc + sparse(rows(:),cols(:),dJ_int_dp(:), (Ci-1)*dim,addim*length(u));
            
            rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim*dim]);
            eqn_p = opt.p_adjt_idx;
            cols  = (Xj_kplus1(s)-1)*xbpdim*addim+repmat((seg.x0_idx(:)-1)'*addim,[dim,pdim]);
            cols  = cols + repelem(eqn_p(:)',dim,dim);
            dJ_int_bc = dJ_int_bc + sparse(rows,cols,-Akplus1sp(:),(Ci-1)*dim,addim*length(u));
        end
            
            dhistdpdp = zeros(dim,pdim,pdim); %==derivative of history at zero
            rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim*pdim]);
            eqn_p = opt.p_adjt_idx;
            cols = (xbpdim*Nj+xbpdim+TDIM)*addim+repmat(eqn_p(:)', dim, pdim)+...
                repelem(0:addim:(pdim-1)*addim, dim,pdim);
            dJ_int_bc = dJ_int_bc + sparse(rows(:),cols(:),dhistdpdp(:),(Ci-1)*dim,addim*length(u));        
        
    else
        
        Xj_kplus1 = Xj_id{k+1};
        
        Sk      = length(Xj_k);
        Skplus1 = length(Xj_kplus1);
        
        X0_adjt_idx = reshape(opt.X0_adjt_idx,[dim,Nj]);
        X1_adjt_idx = reshape(opt.X1_adjt_idx,[dim,Nj]);
        
        if isempty(coupling_seg.coupling{k})
            Ak  = eye(dim);
            Akp = zeros(dim,dim*Sk,pdim);
            Akpp = zeros(dim,dim*Sk,pdim,pdim);
        else
            cfunc = coupling_seg.coupling{k}.func;
            cfunc_dp = coupling_seg.coupling{k}.func_dp;
            cfunc_dpdp = coupling_seg.coupling{k}.func_dpdp;
            
            Mat = cfunc(p);
            Mat_dp = cfunc_dp(p);
            Mat_dpdp = cfunc_dpdp(p);
            rows =  coupling_seg.coupling{k}.rows;
            cols =  coupling_seg.coupling{k}.cols;
            
            if isempty(rows)
            Ak = Mat;
            Akp  = zeros(dim,dim*Sk,pdim);
            Akpp = zeros(dim,dim*Sk,pdim,pdim);    
            else
            Ak  = Mat(rows,cols);   
            Akp  = Mat_dp(rows,cols,1:pdim);   
            Akpp = Mat_dpdp(rows,cols,1:pdim,1:pdim);               
%             Akp  = zeros(dim,dim*Sk,pdim);
%             Akpp = zeros(dim,dim*Sk,pdim,pdim);
            end
            

        end
        
        if isempty(coupling_seg.coupling{k+1})
            Akplus1  = eye(dim);
            Akplus1p = zeros(dim,dim*Skplus1,pdim);
            Akplus1pp = zeros(dim,dim*Skplus1,pdim,pdim);
        else
            cfunc      = coupling_seg.coupling{k+1}.func;
            cfunc_dp   = coupling_seg.coupling{k+1}.func_dp;
            cfunc_dpdp = coupling_seg.coupling{k+1}.func_dpdp;
            rows =  coupling_seg.coupling{k+1}.rows;
            cols =  coupling_seg.coupling{k+1}.cols;
            
            Mat      = cfunc(p);
            Mat_dp   = cfunc_dp(p);
            Mat_dpdp = cfunc_dpdp(p);
            
            if isempty(rows)
             Akplus1   = Mat;   
             Akplus1p  = Mat_dp;
             Akplus1pp = Mat_dpdp;
            else
             Akplus1   = Mat(rows,cols); 
             Akplus1p  = Mat_dp(rows,cols,1:pdim);
             Akplus1pp = Mat_dpdp(rows,cols,1:pdim,1:pdim);
            end
            
   

        end
        
        % contribution w.r.t. x1_jks
        for s=1:Sk
            Aks   = Ak(:,1+(s-1)*dim:s*dim);
            Aksp  = Akp(:,1+(s-1)*dim:s*dim,:);
            Akspp = Akpp(:,1+(s-1)*dim:s*dim,:,:);
            x_jks_bp = reshape(Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim),[dim,NTST*(NCOL+1)]);
            
            % Jacobian w.r.t. p
            rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim*pdim]);
            eqn_p = opt.p_adjt_idx;
            cols = (xbpdim*Nj+xbpdim+TDIM)*addim+repmat(eqn_p(:)', dim, pdim)+...
                repelem(0:addim:(pdim-1)*addim, dim,pdim);
            dJ_int_dp = kron(coupling_seg.id1*x_jks_bp(:),eye(dim))'*reshape(Akspp(:),[dim^2,pdim*pdim]);
            
            dJ_int_bc = dJ_int_bc + sparse(rows(:),cols(:),dJ_int_dp(:), (Ci-1)*dim,addim*length(u));
            
            % Jacobian w.r.t. x1
            rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim*dim]);
            eqn_p = opt.p_adjt_idx;
            cols  = (Xj_k(s)-1)*xbpdim*addim+repmat((seg.x1_idx(:)-1)'*addim,[dim,pdim]);
            cols  = cols + repelem(eqn_p(:)',dim,dim);
            dJ_int_bc = dJ_int_bc + sparse(rows,cols,Aksp(:),(Ci-1)*dim,addim*length(u));
            
        end
        
        % contribution w.r.t. x0_jks
        for s=1:Skplus1
            Akplus1s   = Akplus1(:,1+(s-1)*dim:s*dim);
            Akplus1sp  = Akplus1p(:,1+(s-1)*dim:s*dim,:);
            Akplus1spp = Akplus1pp(:,1+(s-1)*dim:s*dim,:,:);
            x_jkplus1s_bp = reshape(Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim),[dim,NTST*(NCOL+1)]);
            
            % Jacobian w.r.t. p
            rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim*pdim]);
            eqn_p = opt.p_adjt_idx;
            cols = (xbpdim*Nj+xbpdim+TDIM)*addim+repmat(eqn_p(:)', dim, pdim)+...
                repelem(0:addim:(pdim-1)*addim, dim,pdim);
            dJ_int_dp = -kron(coupling_seg.id0*x_jkplus1s_bp(:),eye(dim))'*reshape(Akplus1spp(:),[dim^2,pdim*pdim]);
            
            dJ_int_bc = dJ_int_bc + sparse(rows(:),cols(:),dJ_int_dp(:), (Ci-1)*dim,addim*length(u));
            
            rows = repmat((1+(k-1)*dim:k*dim)', [1,pdim*dim]);
            eqn_p = opt.p_adjt_idx;
            cols  = (Xj_kplus1(s)-1)*xbpdim*addim+repmat((seg.x0_idx(:)-1)'*addim,[dim,pdim]);
            cols  = cols + repelem(eqn_p(:)',dim,dim);
            dJ_int_bc = dJ_int_bc + sparse(rows,cols,-Akplus1sp(:),(Ci-1)*dim,addim*length(u));
        end
        
    end
end

%% Jacobian of Adjoint contributions from mesh conditions on xi_bks and xi_eks
            dJ_xi_ek = []; dJ_xi_bk = [];
            for k=1:Ci-1
                Xj_k     = Xj_id{k};
                if Tj_id{k}~=0
                    Tj_k     = T(Tj_id{k});
                end
                Delta_k  = Delta(Delta_id{k});
                Sk = length(Xj_k);
                
                Ti             = T(Ti_id);
                
                Xj_kplus1      = Xj_id{k+1};
                Tj_kplus1      = T(Tj_id{k+1});
                Delta_kplus1   = Delta(Delta_id{k+1});
                Skplus1        = length(Xj_kplus1);
                
                if Xj_k(1)==0
                    
                    Ti_idx = opt.T_adjt_idx(Ti_id);
                    
                    dTi_dgammaek = 1;
                    dgammaek_dTi = 1;
                    
                    dTi_dgammaekcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltadim)*addim+(2*k-1)*addim+Ti_idx;
                    dgammaek_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + opt.gamma_adjt_idx(2*k);
                    
                    rows = [1;1]; cols = [dTi_dgammaekcols;dgammaek_dTicols];
                    temp = sparse(rows,cols,[dTi_dgammaek;dgammaek_dTi],1,opt.dJcols);
                    dJ_xi_ek = [dJ_xi_ek;temp];
                else
                    gamma_ek = gamma(2*k);
                    
                    dTi_dTjk = -(1/Tj_k^2)*(gamma_ek-Delta_k);
                    dTi_dgammaek = 1/Tj_k;
                    
                    dTi_dDeltak = -(1/Tj_k);
                    
                    dTjk_dTi   = -(1/Tj_k^2)*(gamma_ek-Delta_k);
                    dTjk_dTjk = (2*Ti/Tj_k^3)*(gamma_ek-Delta_k);
                    dTjk_dgammaek = -(Ti/Tj_k^2);
                    
                    dTjk_dDeltak = (Ti/Tj_k^2);
                    
                    dgammaek_dTi   = 1/Tj_k;
                    dgammaek_dTjk = -Ti/Tj_k^2;
                    
                    dDeltak_dTi  = -(1/Tj_k);
                    dDeltak_dTjk =  (Ti/Tj_k^2);
                    
                    eqn_Ti = opt.T_adjt_idx(Ti_id); eqn_Tjk = opt.T_adjt_idx(Tj_id{k}); eqn_gammaek = opt.gamma_adjt_idx(2*k);
                    eqn_Deltak = opt.Delta_adjt_idx(Delta_id{k});
                    
                    dTi_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + eqn_Ti;
                    dTi_dDeltakcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}-1)*addim + eqn_Ti;
                    dTi_dgammaekcols  = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltadim)*addim+(2*k-1)*addim+eqn_Ti;
                    
                    dTjk_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + eqn_Tjk;
                    dTjk_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + eqn_Tjk;
                    dTjk_dDeltakcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k}-1)*addim + eqn_Tjk;
                    dTjk_dgammaekcols = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltadim)*addim+(2*k-1)*addim + eqn_Tjk;
                    
                    dgammaek_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + eqn_gammaek;
                    dgammaek_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + eqn_gammaek;
                    
                    dDeltak_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + eqn_Deltak;
                    dDeltak_dTjkcols = (xbpdim*Nj+xbpdim+Tj_id{k}-1)*addim + eqn_Deltak;
                    
                    rows = ones(11,1);
                    cols = [dTi_dTjkcols;dTi_dDeltakcols;dTi_dgammaekcols;dTjk_dTicols;dTjk_dTjkcols;...
                        dTjk_dDeltakcols; dTjk_dgammaekcols;dDeltak_dTicols;dgammaek_dTicols;dDeltak_dTjkcols;...
                         dgammaek_dTjkcols];
                    vals = [dTi_dTjk;dTi_dDeltak;dTi_dgammaek;dTjk_dTi;dTjk_dTjk;dTjk_dDeltak;...
                        dTjk_dgammaek;dDeltak_dTi;dDeltak_dTjk;dgammaek_dTi;dgammaek_dTjk];
                    temp = sparse(rows,cols,vals,1,opt.dJcols);
                    dJ_xi_ek = [dJ_xi_ek;temp];
                    
                end
                 
                    gamma_bkplus1 = gamma(2*k+1);
                    
                    dTi_dTjkplus1 = -(1/Tj_kplus1^2)*(gamma_bkplus1-Delta_kplus1);
                    dTi_dgammabkplus1 = 1/Tj_kplus1;
                    
                    dTi_dDeltakplus1 = -(1/Tj_kplus1);
                    
                    dTjkplus1_dTi   = -(1/Tj_kplus1^2)*(gamma_bkplus1-Delta_kplus1);
                    dTjkplus1_dTjkplus1 = (2*Ti/Tj_kplus1^3)*(gamma_bkplus1-Delta_kplus1);
                    dTjkplus1_dgammabkplus1 = -(Ti/Tj_kplus1^2);
                    
                    dTjkplus1_dDeltakplus1 = (Ti/Tj_kplus1^2);
                    
                    dgammabkplus1_dTi   = 1/Tj_kplus1;
                    dgammabkplus1_dTjkplus1 = -Ti/Tj_kplus1^2;
                    
                    dDeltakplus1_dTi  = -(1/Tj_kplus1);
                    dDeltakplus1_dTjkplus1 =  (Ti/Tj_kplus1^2);
                    
                    eqn_Ti = opt.T_adjt_idx(Ti_id);
                    eqn_Tjkplus1 = opt.T_adjt_idx(Tj_id{k+1});
                    eqn_Deltakplus1 = opt.Delta_adjt_idx(Delta_id{k+1});
                    eqn_gammabkplus1 = opt.gamma_adjt_idx(2*k+1);
                    
                    
                    dTi_dTjkplus1cols = (xbpdim*Nj+xbpdim+Tj_id{k+1}-1)*addim + eqn_Ti;
                    dTi_dDeltakplus1cols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k+1}-1)*addim + eqn_Ti;
                    dTi_dgammabkplus1cols  = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltadim)*addim+(2*k)*addim+eqn_Ti;
                    
                    dTjkplus1_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + eqn_Tjkplus1;
                    dTjkplus1_dTjkplus1cols = (xbpdim*Nj+xbpdim+Tj_id{k+1}-1)*addim + eqn_Tjkplus1;
                    dTjkplus1_dDeltakplus1cols = (xbpdim*Nj+xbpdim+TDIM+pdim+Delta_id{k+1}-1)*addim + eqn_Tjkplus1;
                    dTjkplus1_dgammbkplus1cols = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltadim)*addim+(2*k)*addim + eqn_Tjkplus1;
                    
                    dgammabkplus1_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + eqn_gammabkplus1;
                    dgammabkplus1_dTjkplus1cols = (xbpdim*Nj+xbpdim+Tj_id{k+1}-1)*addim + eqn_gammabkplus1;
                    
                    dDeltakplus1_dTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + eqn_Deltakplus1;
                    dDeltak_dTjkplus1cols = (xbpdim*Nj+xbpdim+Tj_id{k+1}-1)*addim + eqn_Deltakplus1;
                    
                    rows = ones(11,1);
                    cols = [dTi_dTjkplus1cols;dTi_dDeltakplus1cols;dTi_dgammabkplus1cols;dTjkplus1_dTicols;...
                        dTjkplus1_dTjkplus1cols;dTjkplus1_dDeltakplus1cols;dTjkplus1_dgammbkplus1cols;...
                         dDeltakplus1_dTicols; dDeltak_dTjkplus1cols; dgammabkplus1_dTicols;dgammabkplus1_dTjkplus1cols];
                    vals = [dTi_dTjkplus1;dTi_dDeltakplus1;dTi_dgammabkplus1;dTjkplus1_dTi;...
                        dTjkplus1_dTjkplus1;dTjkplus1_dDeltakplus1;dTjkplus1_dgammabkplus1;...
                        dDeltakplus1_dTi;dDeltakplus1_dTjkplus1; dgammabkplus1_dTi;dgammabkplus1_dTjkplus1];
                    
                    temp = sparse(rows,cols,vals,1,opt.dJcols);
                    dJ_xi_bk = [dJ_xi_bk;temp];
                
                
            end
            


%%
dJ = [dJ_coupling; dJ_fgamma; dJ_int_bc; dJ_xi_ek; dJ_xi_bk];

end

