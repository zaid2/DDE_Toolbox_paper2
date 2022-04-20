function sol = ddaecoll_init_sol(data, t0, x0, y0, p0)
% Copyright (C) Zaid Ahsan, Mingwu Li

t0 = t0(:);
T0 = t0(1);
T  = t0(end)-t0(1); % Interval length
t0 = (t0-t0(1))/T;  % Rescaling - only if T0<>0!
x0 = interp1(t0, x0, data.tbp)'; % Collection of basepoint values

if ~isempty(y0)
y0 = interp1(t0, y0, data.tbp)';  
end

sol.u = [x0(:); y0(:); T0; T; p0(:)];

end
