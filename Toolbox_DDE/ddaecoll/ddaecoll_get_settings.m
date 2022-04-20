function data = ddaecoll_get_settings(prob, tbid, data)
% Copyright (C) Zaid Ahsan, Mingwu Li

defaults.NTST = 10; % Number of mesh intervals
defaults.NCOL = 4;  % Degree of interpolating polynomials
defaults.Dpoints = 'Uniform';
defaults.Dnodes  = 'Gauss';
defaults.Apoints = 'Uniform';
defaults.Anodes  = 'Gauss';
if ~isfield(data, 'ddaecoll')
  data.ddaecoll = [];
end
data.ddaecoll = coco_merge(defaults, coco_merge(data.ddaecoll, ...
  coco_get(prob, tbid))); % Defaults < Stored < User-supplied
if ~coco_exist('TOL', 'class_prop', prob, tbid, '-no-inherit-all')
data.ddaecoll.TOL = coco_get(prob, 'corr', 'TOL')^(2/3);
end
NTST = data.ddaecoll.NTST;
assert(numel(NTST)==1 && isnumeric(NTST) && mod(NTST,1)==0, ...
  '%s: input for option ''NTST'' is not an integer', tbid);
NCOL = data.ddaecoll.NCOL;
assert(numel(NCOL)==1 && isnumeric(NCOL) && mod(NCOL,1)==0, ...
  '%s: input for option ''NCOL'' is not an integer', tbid);

end
