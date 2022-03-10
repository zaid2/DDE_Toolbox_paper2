function [prob, data] = coupling_construct_adjt(prob, tbid, data, sol)
%COLL_CONSTRUCT_ADJT   Add COLL adjoint problem.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coll_construct_opt.m 2872 2015-08-06 20:12:06Z hdankowicz $

seg  = data.coupling_seg.ddaecoll_seg;
opt  = data.coupling_opt;

prob = coco_add_adjt(prob, tbid, @adj, @adj_DU, data, 'aidx', ...
    opt.guidx_adjt,...
    'l0', sol.l0(:), 'tl0', sol.tl0(:), 'adim', opt.adim);

% prob = coco_add_adjt(prob, tbid, @adj, data, 'aidx', ...
%     opt.guidx_adjt,...
%     'l0', sol.l0(:), 'tl0', sol.tl0(:), 'adim', opt.adim);

% prob = coco_add_slot(prob, tbid, @coco_save_data, data, 'save_full');

end

%%
function [data, J] = adj(prob,data,u)

coupling_seg = data.coupling_seg;
seg = coupling_seg.ddaecoll_seg;

NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
dim = seg.dim;
pdim = seg.pdim;
xbpdim = seg.xbpdim;
xcndim = seg.xcndim;

Xibp    = u(coupling_seg.Xibp_idx);
yibp    = u(coupling_seg.yibp_idx);
T       = u(coupling_seg.T_idx);
p       = u(coupling_seg.p_idx);
Deltaj  = u(coupling_seg.Delta_idx);
gamma   = u(coupling_seg.gamma_idx);

tauint = coupling_seg.tauint;

Ci = coupling_seg.Ci;

Xj_id     = coupling_seg.Xj_id;
Tj_id     = coupling_seg.Tj_id;
Ti_id     = coupling_seg.Ti_id;
Deltaj_id = coupling_seg.Deltaj_id;
Nj        = coupling_seg.Nj;
xi        = coupling_seg.xi;

M = seg.M;

%% Adjoint
     
J_xcn   = sparse(xbpdim, Nj*xcndim);
J_T     = zeros(xbpdim,length(coupling_seg.T_idx));
J_p     = zeros(xbpdim,pdim);
J_Delta = zeros(xbpdim, length(coupling_seg.Delta_idx));
J_gamma = sparse(xbpdim, length(coupling_seg.gamma_idx));
 for k=1:Ci
    gammab = gamma(2*k-1);
    gammae = gamma(2*k);
    
    if gammab<0
       gammab=0; 
    end
    
    if gammae>1
       gammae=1; 
    end
    
    idx11   = intersect(find(M>=gammab),find(M<gammae));
    idx12   = intersect(find(M==gammab),find(M==gammae));
    
    idx_ik1 = union(idx11,idx12);
    Mik1    = M(idx_ik1);
    Nik1    = length(Mik1);   
 
    if isempty(Mik1)
        continue;
    else
        
        Xj       = cell2mat(Xj_id(k));
        Sik      = length(Xj);
        Tj       = T(cell2mat(Tj_id(k)));
        Ti       = T(Ti_id);
        Delta    = Deltaj(cell2mat(Deltaj_id(k)));
        
        idx = 1+(idx_ik1(1)-1)*dim:idx_ik1(end)*dim;
        Wik2 = seg.W(idx,:);
        
        Omega_ik = seg.wts2(idx,idx);
        
        if Xj==0
            dhistdt = zeros(dim,Nik1);
            dhistdp = zeros(dim,1:pdim,Nik1);
            
            J_T(:,Ti_id) = J_T(:,Ti_id)-(1/2/NTST)*Wik2'*Omega_ik*diag(repelem(Mik1(:),dim,1))*dhistdt(:);
            J_Delta(:,Deltaj_id{k}) = J_Delta(:,Deltaj_id{k})+(1/2/NTST)*Wik2'*Omega_ik*dhistdt(:);
            J_p(:,1:pdim) = J_p(:,1:pdim)-(1/2/NTST)*Wik2'*Omega_ik*reshape(dhistdp,[Nik1*dim,pdim]);
        else
            if isempty(coupling_seg.coupling{k})
                Aik = eye(dim);
                Bikp = sparse(Nik1*dim, Nik1*dim*Sik*pdim);               
            else               
                cfunc = coupling_seg.coupling{k}.func;
                rows =  coupling_seg.coupling{k}.rows;
                cols =  coupling_seg.coupling{k}.cols;
                Mat = cfunc();
                Aik = Mat(rows,cols);
                
                Aikp = zeros(dim,dim*Sik,pdim);
                temp = repmat(Aikp, [1,Nik1,1]);
                rows = reshape(1:Nik1*dim, [dim, Nik1]);
                rows = repmat(rows, [dim*Sik, pdim]);
                cols = repmat(1:Nik1*dim*Sik*pdim, [dim,1]);
                Bikp = sparse(rows,cols,temp(:));
            end
            
            for s=1:length(Xj)
                %% xjs contribution
                Aiks = Aik(:,1+(s-1)*dim:s*dim);
                
                sb = gammab*Ti/Tj(s)-Delta(s);
                se = gammae*Ti/Tj(s)-Delta(s);
                
                idx_iks = intersect(find(M>=sb),find(M<se));
                Miks1   = M(idx_iks);
                Niks1   = numel(Miks1);
                
                if ~isempty(Miks1)
                    tau_iks1 = (Tj(s)/Ti)*(Miks1(:)+Delta(s));
                    [idx_iks1,tc_iks1] = shift_loc_nodes(tau_iks1,tauint);
                    L1 = coll_L(seg.tm, tc_iks1(:)); pmap1 = kron(L1,eye(dim));
                    L1p = coll_Lp(seg.tm, tc_iks1(:)); dmap1 = kron(L1p,eye(dim));
                    
                    cols = repmat((idx_iks1(:)-1)'*(NCOL+1)*dim, [(NCOL+1)*dim,1])+...
                        repmat((1:(NCOL+1)*dim)', [1,Niks1]);
                    cols = repmat(cols, [dim,1]);
                    rows =  repmat(1:Niks1*dim, (NCOL+1)*dim,1);
                    Wiks1 =  sparse(rows,cols,pmap1',Niks1*dim, NTST*(NCOL+1)*dim);
                    Wiks1p =  sparse(rows,cols,dmap1',Niks1*dim, NTST*(NCOL+1)*dim);
                    
                    idx = 1+(idx_iks(1)-1)*dim+(Xj(s)-1)*xcndim:idx_iks(end)*dim+(Xj(s)-1)*xcndim;
                    J_xcn(:,idx) = J_xcn(:,idx) -(Tj(s)/Ti)*Wiks1'*kron(eye(Niks1),Aiks);
                    
                end
                
                %% Tj_s and Deltaj_s contribution
                tau_iks2 = (Ti/Tj(s))*Mik1(:)-Delta(s);
                [idx_iks2, tc_iks2] = shift_loc_nodes(tau_iks2,tauint);
                
                L2  = coll_L(seg.tm, tc_iks2(:));   pmap2 = kron(L2, eye(dim));
                L2p = coll_Lp(seg.tm, tc_iks2(:));  dmap2 = kron(L2p, eye(dim));
                
                cols = repmat((idx_iks2(:)-1)'*(NCOL+1)*dim, [(NCOL+1)*dim, 1]) +...
                    repmat((1:(NCOL+1)*dim)',[1, Nik1]);
                cols = cols + repmat((Xj(s)-1)*xbpdim,[(NCOL+1)*dim, Nik1]);
                cols = repmat(cols, [dim,1]);
                rows = repmat(1:Nik1*dim, [(NCOL+1)*dim,1]);
                Wiks2  =  sparse(rows,cols,pmap2',Nik1*dim, Nj*xbpdim);
                Wiks2p =  sparse(rows,cols,dmap2',Nik1*dim, Nj*xbpdim);
                
                %       idx = 1+(idx_ik1(1)-1)*dim:idx_ik1(end)*dim;
                J_T(:,Tj_id{k}(s)) = J_T(:,Tj_id{k}(s))+(1/2/NTST)*(Ti/Tj(s)^2)*Wik2'*Omega_ik*...
                    diag(repelem(Mik1(:),dim,1))*kron(eye(Nik1),Aiks)*2*NTST*Wiks2p*Xibp;
                
                J_Delta(:, Deltaj_id{k}(s)) = ...
                    (1/2/NTST)*Wik2'*Omega_ik*kron(eye(Nik1),Aiks)*2*NTST*Wiks2p*Xibp;
            end
            %% Ti contribution
            tau_ik2 = diag(repelem(Mik1(:),Sik,1))*repmat(Ti./Tj(:), [Nik1,1])-...
                repmat(Delta(:),[Nik1,1]);
            [idx_ik2, tc_ik2] = shift_loc_nodes(tau_ik2,tauint);
            L1  = coll_L(seg.tm, tc_ik2(:));   pmap1 = kron(L1, eye(dim));
            L1p = coll_Lp(seg.tm, tc_ik2(:));  dmap1 = kron(L1p, eye(dim));
            cols = repmat((idx_ik2(:)-1)'*(NCOL+1)*dim, [(NCOL+1)*dim, 1]) +...
                repmat((1:(NCOL+1)*dim)',[1, Nik1*Sik]);
            cols = cols + repmat((Xj(:)-1)'*xbpdim,[(NCOL+1)*dim, Nik1]);
            cols = repmat(cols, [dim,1]);
            rows = repmat(1:Nik1*dim*Sik, [(NCOL+1)*dim,1]);
            Wik1  =  sparse(rows,cols,pmap1',Nik1*dim*Sik, Nj*xbpdim);
            Wik1p =  sparse(rows,cols,dmap1',Nik1*dim*Sik, Nj*xbpdim);
            
            G1 = diag(repelem(Mik1(:),Sik,1))*repmat(1./Tj(:),[Nik1,1]);
            G1 = repelem(G1(:), dim, 1);
            
            J_T(:,Ti_id) = J_T(:,Ti_id)-(1/2/NTST)*Wik2'*Omega_ik*kron(eye(Nik1),Aik)*...
                diag(G1(:))*2*NTST*Wik1p*Xibp;
            
            %% p contribution
            J_p(:,1:pdim) = J_p(:,1:pdim) - (1/2/NTST)*Wik2'*Omega_ik*Bikp*kron(eye(pdim),Wik1*Xibp);
            
        end
        

        
    end
 
 end

     rows = 1:(NTST*NCOL+1)*dim; 
     temp = reshape(1:xbpdim,seg.xbp_shp);
     cols = temp(:,seg.tbp_idx);
     
     J_y = sparse(rows,cols,ones((NTST*NCOL+1)*dim,1),(NTST*NCOL+1)*dim,NTST*(NCOL+1)*dim);

     J = [J_xcn J_y' J_T J_p J_Delta J_gamma seg.Q'];

end




%% Jacobian of the Adjoints
function [data,dJ] = adj_DU(prob, data, u)

coupling_seg = data.coupling_seg;
seg = coupling_seg.ddaecoll_seg;
coupling_seg = data.coupling_seg;
opt = data.coupling_opt;

NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
dim = seg.dim;
pdim = seg.pdim;
xbpdim = seg.xbpdim;
xcndim = seg.xcndim;

Xibp    = u(coupling_seg.Xibp_idx);
yibp    = u(coupling_seg.yibp_idx);
T       = u(coupling_seg.T_idx);
p       = u(coupling_seg.p_idx);
Deltaj  = u(coupling_seg.Delta_idx);
gamma   = u(coupling_seg.gamma_idx);

tauint = coupling_seg.tauint;

Ci = coupling_seg.Ci;

Xj_id     = coupling_seg.Xj_id;
Tj_id     = coupling_seg.Tj_id;
Ti_id     = coupling_seg.Ti_id;
Deltaj_id = coupling_seg.Deltaj_id;
Nj        = coupling_seg.Nj;
xi        = coupling_seg.xi;

M = seg.M;

TDIM = length(coupling_seg.T_idx);

addim = opt.addim;
dJ = sparse(xbpdim,addim*length(u));

for k=1:Ci
    gammab = gamma(2*k-1); gammae = gamma(2*k);
    
    idx11  =  intersect(find(M>=gammab),find(M<gammae));
    idx12  =  intersect(find(M==gammab),find(M==gammae));
    idx_ik1 = union(idx11,idx12);
    Mik1   =  M(union(idx11,idx12));
    Nik1   =  length(Mik1);
    
  if isempty(Mik1)
      continue;
  else
    
    Xj       = cell2mat(Xj_id(k));
    Sik      = length(Xj);
    Tj       = T(cell2mat(Tj_id(k)));
    Ti       = T(Ti_id);
    Delta    = Deltaj(cell2mat(Deltaj_id(k)));
    

    
     idx = 1+(idx_ik1(1)-1)*dim:idx_ik1(end)*dim;
     Wik2 = seg.W(idx,:);   
     Omega_ik = seg.wts2(idx,idx);
     
     
     if Xj==0
         dhistdt = zeros(dim,Nik1);
         dhistdp = zeros(dim,1:pdim,Nik1);
         
         dhistdtdt = zeros(dim,Nik1);
         dhistdtdp = zeros(dim,1:pdim,Nik1);
         
            
         dTidTi = -(1/2/NTST)*Wik2'*Omega_ik*diag(repelem(Mik1(:),dim,1))*diag(repelem(Mik1(:),dim,1))*dhistdtdt(:);
         dTidDeltas = (1/2/NTST)*Wik2'*Omega_ik*diag(repelem(Mik1(:),dim,1))*dhistdtdt(:);
         dTidp  = -(1/2/NTST)*Wik2'*Omega_ik*diag(repelem(Mik1(:),dim,1))*reshape(dhistdtdp(:),[Nik1*dim,pdim]);
         
         dDeltasdDeltas = -(1/2/NTST)*Wik2'*Omega_ik*dhistdtdt(:);
         dDeltasdTi = dTidDeltas;
         dDeltasdp  = (1/2/NTST)*Wik2'*Omega_ik*dhistdtdp(:);
         
         dpdTi = dTidp;
         dpdDeltas = dDeltasdp;
%          dpdp = zeros()
         
     else
         
             if isempty(coupling_seg.coupling{k})
        Aik = eye(dim);
        Aikp = zeros(dim,dim*Sik,pdim);
        
%         Bik = sparse(Nik1*dim, Nik1*dim*Sik*pdim);
        Aikpp = zeros(dim,dim*Sik,pdim,pdim);
    else
        
        cfunc = coupling_seg.coupling{k}.func;
        rows =  coupling_seg.coupling{k}.rows;
        cols =  coupling_seg.coupling{k}.cols;
        Mat = cfunc();
        Aik = Mat(rows,cols);
        
        Aikp  = zeros(dim,dim*Sik,pdim);
        Aikpp = zeros(dim,dim*Sik,pdim,pdim);
%         temp = repmat(zeros(dim,Sik*dim,pdim), [1,NTST*NCOL,1]);
%         
%         rows = reshape(1:Nik1*dim,[dim,NTST*NCOL]);
%         rows = repmat(rows,[dim*Sik,pdim]);
%         cols = repmat(1:Nik1*dim*Sik*pdim,[dim,1]);
%         Bik = sparse(rows,cols,temp(:));    
             end
    
         for s=1:length(Xj)
             Aiks  = Aik(:,1+(s-1)*dim:s*dim);
             
             Aiksp = Aikp(:,1+(s-1)*dim:s*dim,:);
             
             %%   Jacobian of Xj_s contribution
             
             sb = gammab*Ti/Tj(s)-Delta(s);
             se = gammae*Ti/Tj(s)-Delta(s);
             
             idx_iks = intersect(find(M>=sb),find(M<se));
             Miks1   = M(idx_iks);
             Niks1   = numel(Miks1);
             
             if ~isempty(Miks1)
                 tau_iks1 = (Tj(s)/Ti)*(Miks1(:)+Delta(s));
                 [idx_iks1,tc_iks1] = shift_loc_nodes(tau_iks1,tauint);
                 L1 = coll_L(seg.tm, tc_iks1(:)); pmap1 = kron(L1,eye(dim));
                 L1p = coll_Lp(seg.tm, tc_iks1(:)); dmap1 = kron(L1p,eye(dim));
                 
                 cols = repmat((idx_iks1(:)-1)'*(NCOL+1)*dim, [(NCOL+1)*dim,1])+...
                     repmat((1:(NCOL+1)*dim)', [1,Niks1]);
                 cols = repmat(cols, [dim,1]);
                 rows =  repmat(1:Niks1*dim, (NCOL+1)*dim,1);
                 Wiks1 =  sparse(rows,cols,pmap1',Niks1*dim, NTST*(NCOL+1)*dim);
                 Wiks1p =  sparse(rows,cols,dmap1',Niks1*dim, NTST*(NCOL+1)*dim);
                 
                 tauplusDelta = diag(repelem(Miks1, dim,1)+Delta(s));
                 dXjsdTi  =  (Tj(s)/Ti^2)*Wiks1'*kron(eye(Niks1),Aiks)-(Tj(s)/Ti)*(-Tj(s)/Ti^2)*2*NTST*Wiks1p'*tauplusDelta*kron(eye(Niks1),Aiks);
                 dXjsdTjs =  (-1/Ti)*Wiks1'*kron(eye(Niks1),Aiks)-(Tj(s)/Ti)*(1/Ti)*2*NTST*Wiks1p'*tauplusDelta*kron(eye(Niks1),Aiks);
                 dXjsdDeltajs =  (-Tj(s)/Ti)*(2*NTST)*(Tj(s)/Ti)*Wiks1p'*kron(eye(Niks1),Aiks);
                 
                 % dXjsdp
                 temp = repmat(Aiksp,[1,Niks1]);
                 rows = reshape(1:Niks1*dim, [dim,Niks1]);
                 rows = repmat(rows, [dim, pdim]);
                 cols = repmat(1:Niks1*dim*pdim, [dim,1]);
                 Biksp = sparse(rows(:),cols(:),temp(:));
                 
                 dXjsdp = -(Tj(s)/Ti)*Wiks1'*Biksp;
                 
                 rows = repmat((1:NTST*(NCOL+1)*dim)',[1,Niks1*dim]);
                 cols = repmat(1+(idx_iks(1)-1)*dim+(Xj(s)-1)*xcndim:idx_iks(end)*dim+(Xj(s)-1)*xcndim,[NTST*(NCOL+1)*dim,1]);
                 
                 dXjsdTjsrows = rows(:);
                 dXjsdTjscols = cols(:) + (xbpdim*Nj+xbpdim+Tj_id{k}(s)-1)*addim;
                 dXjsdTirows = rows(:);
                 dXjsdTicols = cols(:) + (xbpdim*Nj+xbpdim+Ti_id-1)*addim;
                 dXjsdDeltajsrows = rows(:);
                 dXjsdDeltajscols = cols(:) + (xbpdim*Nj+xbpdim+TDIM+pdim+Deltaj_id{k}(s)-1)*addim;
                 
                 dXjsdprows = repmat(rows(:),[1,pdim]);
                 dXjsdpcols = repmat(cols(:),[1,pdim])+(xbpdim*Nj+xbpdim+TDIM)*addim+repmat((0:1:pdim-1)*addim,[xbpdim*Niks1*dim,1]);
                 
                 dXjsrows = [dXjsdTirows(:);dXjsdTjsrows(:);dXjsdprows(:); dXjsdDeltajsrows(:)];
                 dXjscols = [dXjsdTicols(:);dXjsdTjscols(:);dXjsdpcols(:); dXjsdDeltajscols(:)];
                 dXjsvals = [dXjsdTi(:);dXjsdTjs(:);dXjsdp(:);dXjsdDeltajs(:)];
                 
                 dJ = dJ + sparse(dXjsrows,dXjscols,dXjsvals,opt.dJrows,opt.dJcols);
             end
             
             
             
             %%   Jacobian of Tjs,Deltas,Ti contribution
             tau2 = (Ti/Tj(s))*Mik1(:)-Delta(s);
             [idx2, tc2] = shift_loc_nodes(tau2,tauint);
             
             L2   = coll_L(seg.tm, tc2(:));    pmap2 = kron(L2, eye(dim));
             L2p  = coll_Lp(seg.tm, tc2(:));   dmap2 = kron(L2p, eye(dim));
             L2pp = coll_Lpp(seg.tm, tc2(:));  ddmap2 = kron(L2pp, eye(dim));
             
             cols = repmat((idx2(:)-1)'*(NCOL+1)*dim, [(NCOL+1)*dim, 1]) +...
                 repmat((1:(NCOL+1)*dim)',[1, Nik1]);
             cols = cols + repmat((Xj(s)-1)*xbpdim,[(NCOL+1)*dim, Nik1]);
             cols = repmat(cols, [dim,1]);
             rows = repmat(1:Nik1*dim, [(NCOL+1)*dim,1]);
             Wiks2  =  sparse(rows,cols,pmap2',Nik1*dim, Nj*xbpdim);
             Wiks2p =  sparse(rows,cols,dmap2',Nik1*dim, Nj*xbpdim);
             Wiks2pp =  sparse(rows,cols,ddmap2',Nik1*dim, Nj*xbpdim);
             
             %% ===Jacobian of Tjs contribution
             eqn_Tjs = opt.T_adjt_idx(Tj_id{k}(s));
             
             dTjsdXibp = (1/2/NTST)*Wik2'*Omega_ik*(Ti/Tj(s)^2)*...
                 diag(repelem(Mik1(:),dim,1))*kron(eye(Nik1),Aiks)*2*NTST*Wiks2p;
             
             dTjsdXibprows = repmat((1:xbpdim)',[1,Nj*xbpdim]);
             dTjsdXibpcols = repmat(eqn_Tjs,[xbpdim,Nj*xbpdim])+...
                 repmat(0:addim:(Nj*xbpdim-1)*addim,[xbpdim,1]);
             
             H1 = diag(repelem(Mik1(:)/Tj(s),dim,1));
             dTjsdTi = (1/2/NTST)*Wik2'*Omega_ik*(1/Tj(s)^2)*...
                 diag(repelem(Mik1(:),dim,1))*kron(eye(Nik1),Aiks)*2*NTST*Wiks2p*Xibp+...
                 (1/2/NTST)*Wik2'*Omega_ik*(Ti/Tj(s)^2)*...
                 diag(repelem(Mik1(:),dim,1))*kron(eye(Nik1),Aiks)*H1*4*NTST^2*Wiks2pp*Xibp;
             dTjsdTirows = 1:xbpdim;
             dTjsdTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Tjs,[xbpdim,1]);
             
             H2 = diag(repelem(-Ti*Mik1(:)/Tj(s)^2,dim,1));
             dTjsdTjs = (1/2/NTST)*Wik2'*Omega_ik*(-2*Ti/Tj(s)^3)*...
                 diag(repelem(Mik1(:),dim,1))*kron(eye(Nik1),Aiks)*2*NTST*Wiks2p*Xibp+...
                 (1/2/NTST)*Wik2'*Omega_ik*(Ti/Tj(s)^2)*...
                 diag(repelem(Mik1(:),dim,1))*kron(eye(Nik1),Aiks)*H2*4*NTST^2*Wiks2pp*Xibp;
             dTjsdTjsrows = 1:xbpdim;
             dTjsdTjscols = (xbpdim*Nj+xbpdim+Tj_id{k}(s)-1)*addim + repmat(eqn_Tjs,[xbpdim,1]);
             
             dTjsdDeltas = -(1/2/NTST)*Wik2'*Omega_ik*(Ti/Tj(s)^2)*...
                 diag(repelem(Mik1(:),dim,1))*kron(eye(Nik1),Aiks)*4*NTST^2*Wiks2pp*Xibp;
             dTjsdDeltasrows = 1:xbpdim;
             dTjsdDeltascols = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltaj_id{k}(s)-1)*addim + repmat(eqn_Tjs,[xbpdim,1]);
             
             temp = repmat(Aiksp,[1,Nik1]);
             rows = reshape(1:Nik1*dim, [dim,Nik1]);
             rows = repmat(rows, [dim, pdim]);
             cols = repmat(1:Nik1*dim*pdim, [dim,1]);
             Biksp = sparse(rows(:),cols(:),temp(:));
             
             dTjsdp = (1/2/NTST)*(Ti/Tj(s)^2)*Wik2'*Omega_ik*...
                 diag(repelem(Mik1(:),dim,1))*Biksp*kron(eye(pdim),2*NTST*Wiks2p*Xibp);
             dTjsdprows = repmat(1:xbpdim, [1,pdim]);
             dTjsdpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
                 repmat(eqn_Tjs,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
             
             dTjsvals = [dTjsdXibp(:);dTjsdTi(:);dTjsdTjs(:);dTjsdp(:);dTjsdDeltas(:)];
             dTjsrows = [dTjsdXibprows(:);dTjsdTirows(:);dTjsdTjsrows(:);dTjsdprows(:);dTjsdDeltasrows(:)];
             dTjscols = [dTjsdXibpcols(:);dTjsdTicols(:);dTjsdTjscols(:);dTjsdpcols(:);dTjsdDeltascols(:)];
             
             dJ = dJ + sparse(dTjsrows,dTjscols,dTjsvals,opt.dJrows,opt.dJcols);
             
             %% ========Jacobian of Deltas contribution
             eqn_Deltas = opt.Delta_adjt_idx(Deltaj_id{k}(s));
             
             dDeltasdXibp = (1/2/NTST)*Wik2'*Omega_ik*kron(eye(Nik1),Aiks)*2*NTST*Wiks2p;
             rows = 1:xbpdim;
             dDeltasdXibprows = repmat(rows(:),[1,Nj*xbpdim]);
             dDeltasdXibpcols = repmat(eqn_Deltas,[xbpdim,Nj*xbpdim])+...
                 repmat(0:addim:(Nj*xbpdim-1)*addim,[xbpdim,1]);
             
             dDeltasdTjs  = dTjsdDeltas;
             dDeltasdTjsrows = 1:xbpdim;
             dDeltasdTjscols = (xbpdim*Nj+xbpdim+Tj_id{k}(s)-1)*addim + repmat(eqn_Deltas,[xbpdim,1]);
             
             H1 = diag(repelem(Mik1(:)/Tj(s),dim,1));
             dDeltasdTi = (1/2/NTST)*Wik2'*Omega_ik*kron(eye(Nik1),Aiks)*H1*4*NTST^2*Wiks2pp*Xibp;
             dDeltasdTirows = 1:xbpdim;
             dDeltasdTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Deltas,[xbpdim,1]);
             
             dDeltasdDeltas = -(1/2/NTST)*Wik2'*Omega_ik*kron(eye(Nik1),Aiks)*4*NTST^2*Wiks2pp*Xibp;
             dDeltasdDeltasrows = 1:xbpdim;
             dDeltasdDeltascols = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltaj_id{k}(s)-1)*addim + repmat(eqn_Deltas,[xbpdim,1]);
             
             temp = repmat(Aiksp,[1,Nik1]);
             rows = reshape(1:Nik1*dim, [dim,Nik1]);
             rows = repmat(rows, [dim, pdim]);
             cols = repmat(1:Nik1*dim*pdim, [dim,1]);
             Biksp = sparse(rows(:),cols(:),temp(:));
             
             dDeltasdp = (1/2/NTST)*Wik2'*Omega_ik*Biksp*kron(eye(pdim),2*NTST*Wiks2p*Xibp);
             dDeltasdprows = repmat(1:xbpdim, [1,pdim]);
             dDeltasdpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
                 repmat(eqn_Deltas,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
             
             dDeltasvals = [dDeltasdXibp(:);dDeltasdTi(:);dDeltasdTjs(:);dDeltasdp(:);dDeltasdDeltas(:)];
             dDeltasrows = [dDeltasdXibprows(:);dDeltasdTirows(:);dDeltasdTjsrows(:);dDeltasdprows(:);dDeltasdDeltasrows(:)];
             dDeltascols = [dDeltasdXibpcols(:);dDeltasdTicols(:);dDeltasdTjscols(:);dDeltasdpcols(:);dDeltasdDeltascols(:)];
             
             dJ = dJ + sparse(dDeltasrows,dDeltascols,dDeltasvals,opt.dJrows,opt.dJcols);
             
             %% ===Jacobian of Ti contribution
             eqn_Ti = opt.T_adjt_idx(Ti_id);
             
             dTidTjs = dTjsdTi;
             dTidTjsrows = 1:xbpdim;
             dTidTjscols = (xbpdim*Nj+xbpdim+Tj_id{k}(s)-1)*addim + repmat(eqn_Ti,[xbpdim,1]);
             
             dTidDeltas = dDeltasdTi;
             dTidDeltasrows = 1:xbpdim;
             dTidDeltascols = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltaj_id{k}(s)-1)*addim + repmat(eqn_Ti,[xbpdim,1]);
             
             dTidTjsdDeltasrows  = [dTidTjsrows(:);dTidDeltasrows(:)];
             dTidTjsdDeltascols  = [dTidTjscols(:);dTidDeltascols(:)];
             dTidTjsdDeltasvals = [dTidTjs(:);dTidDeltas(:)];
             
             dJ = dJ + sparse(dTidTjsdDeltasrows,dTidTjsdDeltascols,...
                 dTidTjsdDeltasvals,opt.dJrows,opt.dJcols);
             
             %% Jacobian of p contribution
             eqn_p = opt.p_adjt_idx;
             
             dpdTjs  = dTjsdp;
             dpdDeltas = dDeltasdp;
             
             dpdTjsrows = repmat((1:xbpdim)',[1,pdim]);
             dpdTjscols = (xbpdim*Nj+xbpdim+Tj_id{k}(s)-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
             
             dpdDeltasrows = repmat((1:xbpdim)',[1,pdim]);
             dpdDeltascols = (xbpdim*Nj+xbpdim+TDIM+pdim+Deltaj_id{k}(s)-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
             
             dpdTjsdDeltasvals = [dpdTjs(:);dpdDeltas(:)];
             dpdTjsdDeltasrows = [dpdTjsrows(:);dpdDeltasrows(:)];
             dpdTjsdDeltascols = [dpdTjscols(:);dpdDeltascols(:)];
             
             dJ = dJ + sparse(dpdTjsdDeltasrows,dpdTjsdDeltascols,...
                 dpdTjsdDeltasvals,opt.dJrows,opt.dJcols);
             
         end
         
         % dTidXibp dTidTi dTidp
         tau1 = diag(repelem(Mik1(:),Sik,1))*repmat(Ti./Tj(:), [Nik1,1])-...
             repmat(Delta(:),[Nik1,1]);
         [idx1, tc1] = shift_loc_nodes(tau1,tauint);
         
         L1  = coll_L(seg.tm, tc1(:));   pmap1 = kron(L1, eye(dim));
         L1p = coll_Lp(seg.tm, tc1(:));  dmap1 = kron(L1p, eye(dim));
         L1pp = coll_Lpp(seg.tm, tc1(:));  ddmap1 = kron(L1pp, eye(dim));
         
         cols = repmat((idx1(:)-1)'*(NCOL+1)*dim, [(NCOL+1)*dim, 1]) +...
             repmat((1:(NCOL+1)*dim)',[1, Nik1*Sik]);
         cols = cols + repmat((Xj(:)-1)'*xbpdim,[(NCOL+1)*dim, Nik1]);
         cols = repmat(cols, [dim,1]);
         rows = repmat(1:Nik1*dim*Sik, [(NCOL+1)*dim,1]);
         Wik1   =  sparse(rows,cols,pmap1',Nik1*dim*Sik, Nj*xbpdim);
         Wik1p  =  sparse(rows,cols,dmap1',Nik1*dim*Sik, Nj*xbpdim);
         Wik1pp =  sparse(rows,cols,ddmap1',Nik1*dim*Sik, Nj*xbpdim);
         
         G1 = diag(repelem(Mik1(:),Sik,1))*repmat(1./Tj(:),[Nik1,1]);
         G1 = repelem(G1(:), dim, 1);
         dTidXibp = -(1/2/NTST)*Wik2'*Omega_ik*kron(eye(Nik1),Aik)*diag(G1(:))*2*NTST*Wik1p;
         
         rows = 1:xbpdim;
         dTidXibprows = repmat(rows(:),[1,Nj*xbpdim]);
         dTidXibpcols = repmat(eqn_Ti,[xbpdim,Nj*xbpdim])+...
             repmat(0:addim:(Nj*xbpdim-1)*addim,[xbpdim,1]);
         
         dTidTi = -(1/2/NTST)*Wik2'*Omega_ik*kron(eye(Nik1),Aik)*diag(G1(:))*diag(G1(:))*4*NTST^2*Wik1pp*Xibp;
         dTidTirows = 1:xbpdim;
         dTidTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_Ti,[xbpdim,1]);
         
         temp = repmat(Aikp,[1,Nik1]);
         rows = reshape(1:Nik1*dim, [dim,Nik1]);
         rows = repmat(rows, [dim*Sik, pdim]);
         cols = repmat(1:Nik1*dim*Sik*pdim, [dim,1]);
         Bikp = sparse(rows(:),cols(:),temp(:));
         
         dTidp = -(1/2/NTST)*Wik2'*Omega_ik*Bikp*...
             kron(eye(pdim),diag(G1(:))*2*NTST*Wik1p*Xibp);
         dTidprows = repmat(1:xbpdim,[1,pdim]);
         dTidpcols = (xbpdim*Nj+xbpdim+TDIM)*addim +...
             repmat(eqn_Ti,[xbpdim,pdim])+repmat((0:1:pdim-1)*addim,[xbpdim,1]);
         
         dTidXibpdTidprows = [dTidXibprows(:); dTidTirows(:); dTidprows(:)];
         dTidXibpdTidpcols = [dTidXibpcols(:); dTidTicols(:); dTidpcols(:)];
         dTidXibpdTidpvals = [dTidXibp(:); dTidTi(:); dTidp(:)];
         
         dJ = dJ + sparse(dTidXibpdTidprows,dTidXibpdTidpcols,...
             dTidXibpdTidpvals,opt.dJrows,opt.dJcols);
         
         % dpdXibp, dpdTi and dpdp
         dpdXibp = -(1/2/NTST)*Wik2'*Omega_ik*Bikp*kron(eye(pdim),Wik1);
         
         rows = repmat((1:xbpdim)',[1,Nj*xbpdim]);
         dpdXibprows = repmat(rows, [1,pdim]);
         cols = repmat(1:Nj*xbpdim,[xbpdim,1])+repmat(0:addim:(xbpdim*Nj-1)*addim,[xbpdim,1]);
         dpdXibpcols = repmat(cols(:),[1,pdim]) + repmat(eqn_p(:)'-1,[Nj*xbpdim*xbpdim,1]);
         
         dpdTi = dTidp;
         dpdTirows = repmat((1:xbpdim)',[1,pdim]);
         dpdTicols = (xbpdim*Nj+xbpdim+Ti_id-1)*addim + repmat(eqn_p(:)',[xbpdim,1]);
         
         %         dpdp
         temp = repmat(Aikpp, [1,Nik1,1,1]);
         rows = reshape(1:Nik1*dim, [dim,Nik1]);
         rows = repmat(rows, [dim*Sik, pdim^2]);
         cols = repmat(1:dim*Sik*Nik1*pdim^2, [dim,1]);
         Bikpp = sparse(rows(:),cols(:),temp(:));
         
         dpdp  = -(1/2/NTST)*Wik2'*Omega_ik*Bikpp*kron(eye(pdim^2),Wik1*Xibp);
         dpdprows = repmat((1:xbpdim)',[1,pdim*pdim]);
         dpdpcols = (xbpdim*Nj+xbpdim+TDIM)*addim+...
             repmat(eqn_p(:)', xbpdim, pdim)+...
             repelem(0:addim:(pdim-1)*addim, xbpdim,pdim);
         
         dpdXibpdTidp = [dpdXibp(:);dpdTi(:);dpdp(:)];
         dpdXibpdTidprows = [dpdXibprows(:);dpdTirows(:);dpdprows(:)];
         dpdXibpdTidpcols = [dpdXibpcols(:);dpdTicols(:);dpdpcols(:)];
         
         dJ = dJ + sparse(dpdXibpdTidprows,dpdXibpdTidpcols,...
             dpdXibpdTidp,opt.dJrows,opt.dJcols);
     end
     
       

  end
end


% dJ_incorr = reshape(full(dJ), [xbpdim, opt.addim,length(u)]);
% [data, dJ_corr] = coco_ezDFDX('f(o,d,x)', prob, data, @adj, u);
% check = dJ_corr-dJ_incorr;
% max(abs(check(:)))

end









%%
function A=coll_L(ts,tz)

q=numel(ts);
p=numel(tz);

zi=repmat(reshape(tz,[p 1 1]),[1 q q]);
sj=repmat(reshape(ts,[1 q 1]),[p 1 q]);
sk=repmat(reshape(ts,[1 1 q]),[p q 1]);

t1=zi-sk;
t2=sj-sk;
idx=find(abs(t2)<=eps);
t1(idx)=1;
t2(idx)=1;

A = prod(t1./t2,3);
end

function A = coll_Lp(ts, tz)

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1 1]), [1 q q q]);
sj = repmat(reshape(ts, [1 q 1 1]), [p 1 q q]);
sk = repmat(reshape(ts, [1 1 q 1]), [p q 1 q]);
sl = repmat(reshape(ts, [1 1 1 q]), [p q q 1]);

t3 = sj(:,:,:,1)-sk(:,:,:,1);
t4 = zi-sl;
t5 = sj-sl;

idx1 = find(abs(t5)<=eps);
idx2 = find(abs(t3)<=eps);
idx3 = find(abs(sk-sl)<=eps);
t5(union(idx1, idx3)) = 1;
t4(union(idx1, idx3)) = 1;
t3(idx2) = 1;
t3       = 1.0./t3;
t3(idx2) = 0;

A = sum(t3.*prod(t4./t5, 4), 3);

end



function A = coll_Lpp(ts, tz)

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1 1 1]), [1 q q q q]);
sj = repmat(reshape(ts, [1 q 1 1 1]), [p 1 q q q]);
sk = repmat(reshape(ts, [1 1 q 1 1]), [p q 1 q q]);
sl = repmat(reshape(ts, [1 1 1 q 1]), [p q q 1 q]);
sh = repmat(reshape(ts, [1 1 1 1 q]), [p q q q 1]);

t3 = sj(:,:,:,1,1) - sk(:,:,:,1,1);
t4 = sj(:,:,:,:,1) - sl(:,:,:,:,1);
t5 = zi - sh;
t6 = sj - sh;

idx1 = find(abs(t6)<=eps);
idx2 = find(abs(t3)<=eps);
idx3 = find(abs(t4)<=eps);
idx4 = find(abs(sl(:,:,:,:,1)-sk(:,:,:,:,1))<=eps);

idx5 = find(abs(sl-sh)<=eps);
idx6 = find(abs(sk-sh)<=eps);

t3(idx2) = 1;
t3       = 1.0./t3;
t3(idx2) = 0;

t4(union(idx3,idx4)) = 1;
t4                   = 1.0./t4;
t4(union(idx3,idx4)) = 0;

t6(union(idx1, union(idx5, idx6))) = 1;
t5(union(idx1, union(idx5, idx6))) = 1;

A = sum(t3.*sum(t4.*prod(t5./t6,5),4),3);
end

