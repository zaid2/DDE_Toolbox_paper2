function [idx,tc] = shift_loc_nodes(tauh,tauint)
% Copyright (C) Zaid Ahsan
tauh = tauh(:); tauint = tauint(:);
len2 = length(tauint);  %====NTST+1
len1 = length(tauh);    %====Number of shifted time points

taub = repmat(tauh, [1, len2]);
tauintb = repmat(tauint', [len1, 1]);

idx = sum(ceil(heaviside(taub-tauintb)),2);
if idx(1)==0
idx(idx==0)=1;
end
idx(idx==len2)=len2-1;

% NTST = length(tauint)-1;
% idx = floor(interp1(tauint,(1:NTST+1),tauh));

 tc = 2*(len2-1)*(tauh-tauint(idx(:)))-1; 


end
