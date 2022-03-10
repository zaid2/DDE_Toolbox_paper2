function data = msddaebvp_adjt_init_data(prob, tbid, src_data)
%MSBVP_INIT_DATA   Initialize toolbox data for an instance of 'msbvp'.
%
% Populate remaining fields of the toolbox data structure used by 'msbvp'
% function objects.
%
% DATA = MSBVP_INIT_DATA(PROB, TBID, DATA)
%
% DATA - Toolbox data structure.
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: msbvp_init_data.m 2839 2015-03-05 17:09:01Z fschild $

% xnum = 0;
% for i=1:data.nsegs
%   stbid = coco_get_id(tbid, sprintf('seg%d.ddaecoll', i)); % Construct 'coll' toolbox instance identifier
%   fdata = coco_get_func_data(prob, stbid, 'data');     % Extract 'coll' toolbox data
%   xnum  = xnum+numel(fdata.x0_idx);                    % Track total dimension
% end
% 
% data.T0_idx = (1:data.nsegs)';                     % Index array for initial times
% data.T_idx  = data.nsegs+(1:data.nsegs)';          % Index array for interval lengths
% data.x0_idx = 2*data.nsegs+(1:xnum)';              % Index array for trajectory end points at t=0
% data.x1_idx = 2*data.nsegs+xnum +(1:xnum)';        % Index array for trajectory end points at t=1
% data.p_idx  = 2*data.nsegs+2*xnum+(1:fdata.pdim)'; % Index array for problem parameters

data.msbvp_data = src_data;

end
