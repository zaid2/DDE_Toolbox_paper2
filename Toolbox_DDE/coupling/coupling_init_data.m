function data = coupling_init_data(data,p0)
% Copyright (C) Zaid Ahsan

Ci = length(data.Xj);  % number of coupling conditions
data.Ci = Ci;

seg  = data.ddaecoll_seg;
NTST = seg.ddaecoll.NTST; % Number of mesh intervals
NCOL = seg.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = seg.dim;    % State-space dimension
ydim = seg.ydim;
pseg_dim = seg.pdim;      % Number of problem parameters in the segment

pdim = pseg_dim + numel(p0); % Number of prob. parameters in the coupling condition
data.pdim = pdim;
pnew_dim = numel(p0);
data.pnew_dim = pnew_dim;

xbpdim = seg.xbpdim;
ybpdim = seg.ybpdim;

data.tau_pt      = (0:NTST)/NTST;

yiid = data.yiid;     % y indices 
Xj = data.Xj;         % Cell array with each cell containing Sik indices
xi = data.xi;         % Segment for which coupling conditions are written
Tj = data.Tj;

XJ = [];   %===temporary variable: XJ contains all the j indices
TJ = [];   %===temporary varuable: TJ contains all Tjk indices
for k=1:Ci  
    XJ = [XJ; Xj{k}];    
    if Tj{k}~=0
    TJ = [TJ; Tj{k}];
    end
end
Ji = unique(XJ,'sorted');  %======Unique coupling segments 
T_uniq = unique([xi; TJ]);            % unique T values present 

flag = find(Ji==0, 1);     %=====Ji=0 denotes the history segment

if isempty(flag)      % No history segment      
    Nj = length(Ji(:));
    Xibp_idx = repmat(seg.xbp_idx,[1,Nj]) + repmat((Ji(:)-1)'*(xbpdim+ybpdim+2+pseg_dim),[xbpdim,1]); % unique base points present in coupling
    
%     T_uniq = unique([xi; Ji]);            % unique T values present 

elseif flag==1 && length(Ji)==1           % Only history segment present     
    Xibp_idx = []; 
%     T_uniq = xi;
    Ji = -1;        % Ji set to -1 deliberately 
    Nj = 0;         % No xj segments involved
elseif flag==1 && length(Ji)>1            % History segment+Other segments
    Ji = Ji(2:end,1);                     % Ignoring the history
    Nj = length(Ji(:)); 
    data.Nj = Nj;
    Xibp_idx = repmat(seg.xbp_idx,[1,Nj]) + repmat((Ji-1)'*(xbpdim+ybpdim+2+pseg_dim),[xbpdim,1]);
    
%     T_uniq = unique([xi; Ji]);
end

data.T_uniq = T_uniq;  
               

% Identifiers to extract indices from the unique arrays
%==1. Xj_id{k}: Contains the identifier to the xjs segments
%==2. Tj_id{k}: Contains the identifier to extract Tj_ks from the array T
%==3. Ti_id:    Contains the identifier to extract Ti from the array T
%==4. Delta_id{k}: Contains the identifier to extract Delta_ks from the array Delta

[~,Ti_id] = ismember(xi, T_uniq);     % Ti_id: location of Ti in the array T_uniq 

for k=1:Ci
    [~,Xj_id{k}] = ismember(Xj{k}, Ji); % Special case when Ji=-1, Xj_id{k}=0. This helps to identify the presence of history segment during collocation
    [~,Tj_id{k}] = ismember(Tj{k}, T_uniq);
    if k==1
        Delta_id{k}  = 1;
    else
        Delta_id{k}  = Delta_id{k-1}(end) + 1;
    end
end

% Identifiers below are used during imposition of the coupling conditions
data.Ji = Ji;
data.Xj_id = Xj_id;
data.Tj_id = Tj_id;
data.Ti_id = Ti_id;
data.Deltaj_id = Delta_id;

temp = repmat(yiid(:),[1,NTST*(NCOL+1)])+repmat(0:ydim:ydim*(NTST*(NCOL+1)-1), [dim,1]); % Multiple delays
yibp_idx = seg.ybp_idx(temp(:)) + (xi-1)*(xbpdim+ybpdim+2+pseg_dim);
Tidx   = seg.T_idx + (T_uniq(:)-1)*(xbpdim+ybpdim+2+pseg_dim);
pidx   = seg.p_idx;

if ~isempty(Xibp_idx)
guidx = [Xibp_idx(:); yibp_idx(:); Tidx(:); pidx(:)];
else
guidx = [yibp_idx(:); Tidx(:); pidx(:)];    % Only history segment is present
end

data.guidx = guidx;


if ~isempty(Xibp_idx)
data.Xibp_idx  = (1:Nj*xbpdim)';   
data.yibp_idx  = (data.Xibp_idx(end)+1:data.Xibp_idx(end)+xbpdim)';
data.Nj = Nj;
else
data.Xibp_idx  = [];
data.yibp_idx  = (1:xbpdim)';
data.Nj = 0;
end

data.T_idx     = data.yibp_idx(end)+(1:length(Tidx))';
data.p_idx     = data.T_idx(end) + (1:(pseg_dim+pnew_dim))';
data.Delta_idx = data.p_idx(end) + (1:Ci)';
data.gamma_idx = data.Delta_idx(end) + (1:2*Ci)';

if ~isempty(p0)
data.pnew_idx  = data.T_idx(end) + pseg_dim + (1:pnew_dim)';
else
data.pnew_idx = [];    
end

% data.p_idx = [data.pseg_idx(:);data.pnew_idx(:)];

data.id0 = [eye(dim); zeros(xbpdim-dim,dim)]';
data.id1 = [zeros(xbpdim-dim,dim); eye(dim)]';

% dim_xieks = 0; dim_xibks = 0;
% for k=1:Ci-1
%    Xj_k = Xj_id{k};
%    Sk = length(Xj_k);
%    dim_xieks = dim_xieks + Sk;
%    
%    Xj_kplus1 = Xj_id{k+1};
%    Skplus1 = length(Xj_kplus1);
%    dim_xibks = dim_xibks + Skplus1;
% end

data.dim_xiek = Ci-1;
data.dim_xibk = Ci-1;


end






