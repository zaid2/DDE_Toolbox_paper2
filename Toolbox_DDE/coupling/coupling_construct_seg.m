function prob = coupling_construct_seg(prob, tbid, data,sol)
% Copyright (C) Zaid Ahsan

prob = coco_add_func(prob, tbid, @coup_dde_F, @coup_dde_DFDU, data,...
                                     'zero', 'uidx', data.guidx, 'u0',sol);

if ~isempty(data.pnames) % Optional monitor functions
  fid  = coco_get_id(tbid, 'pars');
  uidx = coco_get_func_data(prob, tbid, 'uidx'); % Context-dependent index array
  prob = coco_add_pars(prob, fid, uidx(data.pnew_idx), data.pnames);
end

prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end


function [data,f] = coup_dde_F(prob, data, u)

seg = data.ddaecoll_seg;

NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
dim = seg.dim;
xbpdim = seg.xbpdim;
xcndim = seg.xcndim;

if ~isempty(data.Xibp_idx)
Xibp = u(data.Xibp_idx);
else
Xibp = [];        
end

ybp = u(data.yibp_idx);
T    = u(data.T_idx);
p    = u(data.p_idx);
Delta = u(data.Delta_idx);
gamma = u(data.gamma_idx);

tbp = seg.tbp;  % without duplication

tau_pt = data.tau_pt;
Ci = data.Ci;

Nj   = data.Nj;   % total number of unique j's for the ith segment

Xj_id = data.Xj_id;
Tj_id = data.Tj_id;
Ti_id = data.Ti_id;
Delta_id = data.Deltaj_id;

%% Coupling conditions
f_coupling = zeros(NTST*(NCOL+1)*dim,1);
bp = 0; 


kbp = floor(interp1([gamma(1:2:end-1);1],1:Ci+1,tbp));
kbp(kbp>Ci) = Ci;

for k = 1:Ci
    Xj_k     = Xj_id{k};
    Sk       = length(Xj_k);
    Ti       = T(Ti_id);
    if Tj_id{k}~=0
    Tj_k     = T(Tj_id{k});
    end
    Delta_k  = Delta(Delta_id{k});
        
    idx_k   = find(kbp==k);
    tbp_k   = tbp(idx_k);
    Nbp_k   = numel(tbp_k);
    
    
    if isempty(tbp_k)
        continue;
    else
        
        if Xj_k(1)==0
            history = data.coupling{k}.func;
            t = Ti*tbp_k-Delta_k(1);
            fhist = history(t(:)',repmat(p,[1,Nbp_k]));
            
            temp = reshape(1:xbpdim,seg.xbp_shp); %===for y variable to get (Nm+1)n index
            idx = temp(:,idx_k);
            ybp_k = ybp(idx(:)); 
            f_coupling(bp+1:bp+Nbp_k*dim,1) = ybp_k-fhist(:);
            
        else
            if isempty(data.coupling{k})
                Ak  = eye(dim);
            else
                cfunc = data.coupling{k}.func;
                Mat = cfunc(p);
                rows =  data.coupling{k}.rows;
                cols = data.coupling{k}.cols;
                
                if isempty(rows)
                    Ak = Mat;
                else
                    Ak = Mat(rows,cols);
                end
                
            end
            
            taubp_lshift = (Ti/Tj_k)*(tbp_k(:)-kron(ones(Nbp_k,1),Delta_k));
            mesh = 0:1/NTST:1;
            jbp_lshift = floor(interp1(mesh,1:NTST+1,taubp_lshift,'linear','extrap'));
            jbp_lshift(jbp_lshift>NTST) = NTST;
            jbp_lshift(jbp_lshift<=0)    = 1;
            
            rows = repmat((jbp_lshift(:)-1)'*(NCOL+1), [(NCOL+1),1]) +...
                repmat((1:(NCOL+1))',[1, Nbp_k]);
            cols = repmat(1:Nbp_k,[NCOL+1,1]);
            
            tc = 2*NTST*taubp_lshift(:,1)+1-2*jbp_lshift(:,1);
            Lc = coll_L(seg.tm, tc);
            Lbp_lshift = sparse(rows,cols,Lc',(NCOL+1)*NTST,Nbp_k);
            
            Sx = 0;    %==temporary variable to add contributions from different s
                
            for s=1:Sk
                Aks = Ak(:,1+(s-1)*dim:s*dim);
                x_jks_bp = Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim);
                Sx = Sx + kron(Lbp_lshift',Aks)*x_jks_bp(:);
            end
            
            temp = reshape(1:xbpdim,seg.xbp_shp); %===for y variable to get (Nm+1)n index
            idx = temp(:,idx_k);
            ybp_k = ybp(idx(:));   
            f_coupling(bp+1:bp+Nbp_k*dim,1) = ybp_k-Sx;
        end
    end
    bp = bp+Nbp_k*dim;
end

%% Mesh conditions on gamma
f_gamma = zeros(Ci+1,1);
for k=1:Ci+1
    if k==1
       gamma_b1 = gamma(1); 
       f_gamma(k,1) = gamma_b1;
    elseif 1<k && k<Ci+1
       gamma_ekminus1 = gamma(2*k-2); gamma_bk = gamma(2*k-1); 
       f_gamma(k,1) = gamma_ekminus1-gamma_bk;
    elseif k==Ci+1
       gamma_eCi = gamma(end); 
       f_gamma(k,1)  = gamma_eCi-1; 
    end
end

%% ====Interior point boundary conditions
f_int_bc = zeros((Ci-1)*dim,1);
for k=1:Ci-1
  %%  
     Xj_k = Xj_id{k};
     
        if Xj_k(1)==0
            history = data.coupling{k}.func;
            if isempty(data.coupling{k+1})
                Akplus1  = eye(dim);
            else
                cfunc = data.coupling{k+1}.func;
                rows =  data.coupling{k+1}.rows;
                cols = data.coupling{k+1}.cols;
                
                Mat = cfunc(p);
                
                if isempty(rows)
                Akplus1 = Mat;
                else
                Akplus1 = Mat(rows,cols);    
                end
                
            end
            
            Xj_kplus1 = Xj_id{k+1};
            Skplus1 = length(Xj_kplus1);
            Sx0 = 0;
            for s=1:Skplus1
                Akplus1s = Akplus1(:,1+(s-1)*dim:s*dim);
                x_jkplus1s_bp = Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim);
                Sx0 = Sx0 + Akplus1s*data.id0*x_jkplus1s_bp(:);
            end          
            f_int_bc(1+(k-1)*dim:k*dim,1) =  history(0,p)-Sx0;        
            
        else
  %%          
        if isempty(data.coupling{k})
            Ak  = eye(dim);
        else
            cfunc = data.coupling{k}.func;
            rows =  data.coupling{k}.rows;
            cols = data.coupling{k}.cols;
            Mat = cfunc(p);
            
            if isempty(rows)
                Ak = Mat;
            else
                Ak = Mat(rows,cols);
            end
            
            
        end
        
        if isempty(data.coupling{k+1})
            Akplus1  = eye(dim);
        else
            cfunc = data.coupling{k+1}.func;
            rows =  data.coupling{k+1}.rows;
            cols = data.coupling{k+1}.cols;
            Mat = cfunc(p);
            
            if isempty(rows)
            Akplus1 = Mat;
            else
            Akplus1 = Mat(rows,cols);            
            end
        end
        
        Xj_k = Xj_id{k};
        Sk = length(Xj_k);
        
        Sx1 = 0;
        for s=1:Sk
        Aks = Ak(:,1+(s-1)*dim:s*dim);
        x_jks_bp = Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim);
        Sx1 = Sx1 + Aks*data.id1*x_jks_bp(:);
        end
        
        Xj_kplus1 = Xj_id{k+1};
        Skplus1 = length(Xj_kplus1);
        Sx0 = 0;
        for s=1:Skplus1
        Akplus1s = Akplus1(:,1+(s-1)*dim:s*dim);
        x_jkplus1s_bp = Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim);
        Sx0 = Sx0 + Akplus1s*data.id0*x_jkplus1s_bp(:);
        end
        
        f_int_bc(1+(k-1)*dim:k*dim,1) =  Sx1-Sx0;
        end 
end    

%% Mesh conditions on xi_biks and xi_eiks
f_xi_ek = [];
f_xi_bk = [];
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
        gamma_ek = gamma(2*k);
        xi_eks = Ti*gamma_ek-Delta_k(s);
        f_xi_ek = [f_xi_ek; xi_eks];
    else  
            gamma_ek = gamma(2*k);
            xi_ek = (Ti/Tj_k)*(gamma_ek-Delta_k)-1;
            f_xi_ek = [f_xi_ek; xi_ek];

    end
            gamma_bkplus1 = gamma(2*k+1);
            xi_bkplus1   = (Ti/Tj_kplus1)*(gamma_bkplus1-Delta_kplus1);
            f_xi_bk = [f_xi_bk; xi_bkplus1];
        
end


%%      Zero problem

        f = [f_coupling; f_gamma; f_int_bc; f_xi_ek; f_xi_bk];
    
end



function [data,J] = coup_dde_DFDU(prob,data,u)
seg = data.ddaecoll_seg;

NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
dim = seg.dim;
xbpdim = seg.xbpdim;
xcndim = seg.xcndim;

if ~isempty(data.Xibp_idx)
Xibp = u(data.Xibp_idx);
Xbpdim    = length(Xibp);
else
Xibp = [];        
end

ybp = u(data.yibp_idx);
T    = u(data.T_idx);
p    = u(data.p_idx);
Delta = u(data.Delta_idx);
gamma = u(data.gamma_idx);

tbp = seg.tbp;  % without duplication

tau_pt = data.tau_pt;
Ci = data.Ci;

Nj   = data.Nj;   % total number of unique j's for the ith segment

Xj_id = data.Xj_id;
Tj_id = data.Tj_id;
Ti_id = data.Ti_id;
Delta_id = data.Deltaj_id;

cntdim = dim*(NTST-1); 


Tdim      = length(T);
Deltadim  = length(Delta);
Gammadim  = length(gamma);
pdim      = length(p);
% pseg_dim  = seg.pdim;
% pnew_dim =  data.pnew_dim;

%% Jacobian of coupling condition 
if ~isempty(data.Xibp_idx)
    dfcoupling_dx = sparse(NTST*(NCOL+1)*dim,Xbpdim);
else
    dfcoupling_dx = [];
end
dfcoupling_dT        = zeros(NTST*(NCOL+1)*dim, Tdim);
dfcoupling_dDelta    = zeros(NTST*(NCOL+1)*dim, Deltadim);
dfcoupling_dp        = zeros(NTST*(NCOL+1)*dim, pdim);

bp = 0; 

kbp = floor(interp1([gamma(1:2:end-1);1],1:Ci+1,tbp));
kbp(kbp>Ci) = Ci;

for k=1:Ci
    Xj_k    = Xj_id{k};
    Sk      = length(Xj_k);
    Ti      = T(Ti_id);
    if Tj_id{k}~=0
    Tj_k    = T(Tj_id{k});
    end
    Delta_k = Delta(Delta_id{k});
    
    idx_k = find(kbp==k);
    tbp_k  = tbp(idx_k);
    Nbp_k  = numel(tbp_k);
    
    if isempty(tbp_k)
        continue;
    else
        %%
        if Xj_k(1)==0
%             dhistdt = zeros(dim,Nbp_k);
%             dhistdp = zeros(dim,pdim,Nbp_k);
            t = Ti*tbp_k-Delta_k(1);
            dhistdt = data.coupling{k}.func_dt(t(:)',repmat(p,[1,Nbp_k]));
            dhistdp = data.coupling{k}.func_dp(t(:)',repmat(p,[1,Nbp_k]));
            
            dfcoupling_dT(bp+1:bp+Nbp_k*dim, Ti_id) = -kron(diag(tbp_k(:)),eye(dim))*dhistdt(:);
            dfcoupling_dDelta(bp+1:bp+Nbp_k*dim, Delta_id{k}) = dhistdt(:);
            
            dprows = repmat(reshape(1:Nbp_k*dim, [dim Nbp_k]), [pdim 1]); % Index array for vectorization
            dpcols = repmat(1:pdim, [dim Nbp_k]);                         % Index array for vectorization
            dhistdp = sparse(dprows(:),dpcols(:),dhistdp(:));
            
            dfcoupling_dp(bp+1:bp+Nbp_k*dim,:) = -dhistdp(:,1:pdim);

        else
            %%
            
            if isempty(data.coupling{k})
                Ak  = eye(dim);
                Akp = zeros(dim,dim*Sk,pdim);
            else
                cfunc    = data.coupling{k}.func;
                cfunc_dp = data.coupling{k}.func_dp;
                Mat    = cfunc(p);
                Mat_dp = cfunc_dp(p);
                
                rows =  data.coupling{k}.rows;
                cols = data.coupling{k}.cols;
                
                if isempty(rows)
                Ak  = Mat;
                Akp = Mat_dp;
                else
                Ak  = Mat(rows,cols);
                Akp = Mat_dp(rows,cols,1:pdim);
                end
                  
            end
            
            taubp_lshift = (Ti/Tj_k)*(tbp_k(:)-kron(ones(Nbp_k,1),Delta_k));
            mesh = 0:1/NTST:1;
            jbp_lshift = floor(interp1(mesh,1:NTST+1,taubp_lshift,'linear','extrap'));
            jbp_lshift(jbp_lshift>NTST) = NTST;
            jbp_lshift(jbp_lshift<0) =  1;
            
            rows = repmat((jbp_lshift(:)-1)'*(NCOL+1), [(NCOL+1),1]) +...
                repmat((1:(NCOL+1))',[1, Nbp_k]);
            cols = repmat(1:Nbp_k,[NCOL+1,1]);
            tc = 2*NTST*taubp_lshift(:,1)+1-2*jbp_lshift(:,1);
            Lc = coll_L(seg.tm, tc); Lcp = coll_Lp(seg.tm, tc);
            Lbp_lshift  = sparse(rows,cols,Lc',(NCOL+1)*NTST,Nbp_k);
            Lbp_lshiftp = sparse(rows,cols,Lcp',(NCOL+1)*NTST,Nbp_k);
                
            for s=1:Sk
                Aks = Ak(:,1+(s-1)*dim:s*dim);
                Aksp = Akp(:,1+(s-1)*dim:s*dim,1:pdim);
                
                idx = 1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim;
                dfcoupling_dx(bp+1:bp+Nbp_k*dim,idx) = - kron(Lbp_lshift',Aks);
                
                x_jks_bp = reshape(Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim),[dim,NTST*(NCOL+1)]);
                
                dfcoupling_dT(bp+1:bp+Nbp_k*dim, Ti_id) =...
                    dfcoupling_dT(bp+1:bp+Nbp_k*dim, Ti_id)-2*NTST*...
                    (1/Tj_k)*kron((Lbp_lshiftp*diag(tbp_k(:)-Delta_k))',Aks)*x_jks_bp(:);
                
                dfcoupling_dT(bp+1:bp+Nbp_k*dim, Tj_id{k}) =...
                    dfcoupling_dT(bp+1:bp+Nbp_k*dim, Tj_id{k})-2*NTST*...
                    (-Ti/Tj_k^2)*kron((Lbp_lshiftp*diag(tbp_k(:)-Delta_k))',Aks)*x_jks_bp(:);
                
                dfcoupling_dDelta(bp+1:bp+Nbp_k*dim, Delta_id{k}) = ...
                     dfcoupling_dDelta(bp+1:bp+Nbp_k*dim, Delta_id{k})-...
                                       (-Ti/Tj_k)*2*NTST*kron(Lbp_lshiftp',Aks)*x_jks_bp(:);
                
                                   
                                   
                for i=1:pdim
                    dfcoupling_dp(bp+1:bp+Nbp_k*dim,i) = -kron(Lbp_lshift',Aksp(:,:,i))*x_jks_bp(:) + dfcoupling_dp(bp+1:bp+Nbp_k*dim,i);
                end
                
            end
        end
    end
    bp = bp+Nbp_k*dim;
end


J_coupling = [dfcoupling_dx eye(xbpdim) dfcoupling_dT dfcoupling_dp dfcoupling_dDelta sparse(NTST*(NCOL+1)*dim,Gammadim)];



%% Jacobian of continuity conditions for ybp


%% Jacobian of Mesh conditions on gamma
J_fgamma = sparse(Ci+1,length(u));
gamma_idx = data.gamma_idx;
for k=1:Ci+1
    if k==1
        J_fgamma(1,gamma_idx(1)) = 1;
    elseif 1<k && k<Ci+1
        J_fgamma(k,gamma_idx(2*k-2)) =  1;
        J_fgamma(k,gamma_idx(2*k-1)) = -1;
    elseif k==Ci+1
        J_fgamma(k,gamma_idx(end)) = 1;
    end
end


%% Jacobian of internal boundary conditions
J_int_bc = sparse((Ci-1)*dim,length(u));
for k=1:Ci-1
 Xj_k = Xj_id{k};
 Xj_kplus1 = Xj_id{k+1};
 
    if Xj_k(1)==0
%%               
        Xj_kplus1 = Xj_id{k+1};
        Skplus1 = length(Xj_kplus1);
        if isempty(data.coupling{k+1})
            Akplus1  = eye(dim);
            Akplus1p = zeros(dim,Skplus1*dim,pdim);
        else
            cfunc    = data.coupling{k+1}.func;
            cfunc_dp = data.coupling{k+1}.func_dp;
            
            Mat    = cfunc(p);
            Mat_dp = cfunc_dp(p);
            
            rows =  data.coupling{k+1}.rows;
            cols = data.coupling{k+1}.cols;
            
            if isempty(rows)
              Akplus1 = Mat;  
              Akplus1p = Mat_dp;
            else
              Akplus1  = Mat(rows,cols);   
              Akplus1p = Mat_dp(rows,cols,1:pdim);  
            end
        end
        
        for s=1:Skplus1
            Akplus1s  = Akplus1(:,1+(s-1)*dim:s*dim);
            Akplus1sp = Akplus1p(:,1+(s-1)*dim:s*dim,:);
            
            % w.r.t. base points
            idx = 1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim;
            J_int_bc(1+(k-1)*dim:k*dim,idx)= J_int_bc(1+(k-1)*dim:k*dim,idx)-Akplus1s*data.id0;
            
            x_jkplus1s_bp = Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim);
            J_int_dp = -kron(data.id0*x_jkplus1s_bp,eye(dim))'*reshape(Akplus1sp(:),[dim^2,pdim]);
            
            J_int_bc(1+(k-1)*dim:k*dim,data.p_idx) = J_int_bc(1+(k-1)*dim:k*dim,data.p_idx) + J_int_dp;
        end
        
%         dhistdp = zeros(dim,pdim);
        dhistdp = data.coupling{k}.func_dp(0,p);
        J_int_bc(1+(k-1)*dim:k*dim,data.p_idx) = J_int_bc(1+(k-1)*dim:k*dim,data.p_idx)+dhistdp(:,1:pdim);
        
      
    else
 %%      
        Sk = length(Xj_k);
        Skplus1 = length(Xj_kplus1);
        
        if isempty(data.coupling{k})
            Ak  = eye(dim);
            Akp = zeros(dim,Sk*dim,pdim);
        else
            cfunc    = data.coupling{k}.func;
            cfunc_dp = data.coupling{k}.func_dp;
            Mat    = cfunc(p);
            Mat_dp = cfunc_dp(p);
            
            rows =  data.coupling{k}.rows;
            cols = data.coupling{k}.cols;
            
            if isempty(rows)
                Ak  = Mat;
                Akp = Mat_dp;
            else
                Ak  = Mat(rows,cols);
                Akp = Mat_dp(rows,cols,1:pdim);
            end           
        end
        
        if isempty(data.coupling{k+1})
            Akplus1  = eye(dim);
            Akplus1p = zeros(dim,Skplus1*dim,pdim);
        else
            cfunc    = data.coupling{k+1}.func;
            cfunc_dp = data.coupling{k+1}.func_dp;
            
            Mat = cfunc(p);
            Mat_dp = cfunc_dp(p);
            
            rows =  data.coupling{k+1}.rows;
            cols = data.coupling{k+1}.cols;
            
            if isempty(rows)
               Akplus1  = Mat;
               Akplus1p = Mat_dp;
            else
               Akplus1  = Mat(rows,cols);
               Akplus1p = Mat_dp(rows,cols,1:pdim);
%                Akplus1p = zeros(dim,Skplus1*dim,pdim);
            end
            
%             Mat = cfunc(p);
%             Akplus1 = Mat(rows,cols);
%             Akplus1p = zeros(dim,Skplus1*dim,pdim);
        end
        
      
        
        for s=1:Sk
            Aks = Ak(:,1+(s-1)*dim:s*dim);
            Aksp = Akp(:,1+(s-1)*dim:s*dim,:);
            idx = 1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim;
            J_int_bc(1+(k-1)*dim:k*dim,idx)= J_int_bc(1+(k-1)*dim:k*dim,idx)+Aks*data.id1;
            
            x_jks_bp = Xibp(1+(Xj_k(s)-1)*xbpdim:Xj_k(s)*xbpdim);
            J_int_dp = kron(data.id1*x_jks_bp,eye(dim))'*reshape(Aksp(:),[dim^2,pdim]);
            
            J_int_bc(1+(k-1)*dim:k*dim,data.p_idx) = J_int_bc(1+(k-1)*dim:k*dim,data.p_idx) + J_int_dp(:,1:pdim);
            
        end
        
        Xj_kplus1 = Xj_id{k+1};
        Skplus1 = length(Xj_kplus1);
        for s=1:Skplus1
            Akplus1s = Akplus1(:,1+(s-1)*dim:s*dim);
            Akplus1sp = Akplus1p(:,1+(s-1)*dim:s*dim,:);
            idx = 1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim;
            J_int_bc(1+(k-1)*dim:k*dim,idx)= J_int_bc(1+(k-1)*dim:k*dim,idx)-Akplus1s*data.id0;
            
            x_jkplus1s_bp = Xibp(1+(Xj_kplus1(s)-1)*xbpdim:Xj_kplus1(s)*xbpdim);
            J_int_dp = -kron(data.id0*x_jkplus1s_bp,eye(dim))'*reshape(Akplus1sp(:),[dim^2,pdim]);
            
            J_int_bc(1+(k-1)*dim:k*dim,data.p_idx) = J_int_bc(1+(k-1)*dim:k*dim,data.p_idx) + J_int_dp(:,1:pdim);
            
        end
    end
end



%% Mesh conditions on xi_bks and xi_eks
J_xi_ek = [];
J_xi_bk = [];
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
        gamma_ek = gamma(2*k);
%         xi_eiks = Ti*gamma_eik-Delta_k(s);    
        dxiek_dTi = gamma_ek;
        dxiek_dgammaek = Ti;
        dxiek_dDeltak  = -1;
        
        Ti_idx = data.T_idx(Ti_id);
        gammaek_idx = data.gamma_idx(2*k);
        Deltak_idx  = data.Delta_idx(Delta_id{k}(1));
        
        rows = ones(3,1); 
        cols = [Ti_idx;gammaek_idx;Deltak_idx];
        vals = [dxiek_dTi; dxiek_dgammaek; dxiek_dDeltak];
        temp = sparse(rows,cols,vals,1,length(u));
        J_xi_ek = [J_xi_ek;temp];
        else
%         xi_eiks = (Ti/Tj_k(s))*gamma_eik-Delta_k(s);
        gamma_ek = gamma(2*k);
        dxiek_dTi = (1/Tj_k)*(gamma_ek-Delta_k);
        dxieik_dTjk = -(Ti/Tj_k^2)*(gamma_ek-Delta_k);
        dxiek_dgammaek = Ti/Tj_k;
        dxiek_dDeltak  = -(Ti/Tj_k);
        
        Ti_idx = data.T_idx(Ti_id);
        Tjk_idx = data.T_idx(Tj_id{k});
        gammaek_idx = data.gamma_idx(2*k);
        Deltak_idx  = data.Delta_idx(Delta_id{k});
        
        rows = ones(4,1);
        cols = [Ti_idx;Tjk_idx; gammaek_idx;Deltak_idx];
        vals = [dxiek_dTi dxieik_dTjk dxiek_dgammaek dxiek_dDeltak];
        temp = sparse(rows,cols,vals,1,length(u));
        
        J_xi_ek = [J_xi_ek;temp];
        
    end
    
    gamma_bkplus1 = gamma(2*k+1);
    
    dxibkplus1_dTi            = (1/Tj_kplus1)*(gamma_bkplus1-Delta_kplus1);
    dxibkplus1_dTjkplus1      = -(Ti/Tj_kplus1^2)*(gamma_bkplus1-Delta_kplus1);
    dxibkplus1_dgammabkplus1  = (Ti/Tj_kplus1);
    dxibkplus1_dDeltakplus1   = -(Ti/Tj_kplus1);
    
    Ti_idx = data.T_idx(Ti_id);
    Tjk_idx = data.T_idx(Tj_id{k+1});
    gammabkplus1_idx = data.gamma_idx(2*k+1);
    Deltakplus1_idx  = data.Delta_idx(Delta_id{k+1});
    
    rows = ones(4,1);
    cols = [Ti_idx;Tjk_idx; gammabkplus1_idx;Deltakplus1_idx];
    vals = [dxibkplus1_dTi dxibkplus1_dTjkplus1 dxibkplus1_dgammabkplus1 dxibkplus1_dDeltakplus1];
    temp = sparse(rows, cols,vals,1,length(u));
    J_xi_bk = [J_xi_bk;temp];
        
end


%% Overall Jacobian

J = [J_coupling;J_fgamma;J_int_bc;J_xi_ek;J_xi_bk];


end
































