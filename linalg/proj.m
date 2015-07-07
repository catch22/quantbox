function rho = proj(indices, dims)
% Return rank-one projection matrix onto computational basis vector |indices>.
%
% Usage
% =====
%
% RHO = proj(INDICES, DIMS)
%
%
% Examples
% ========
%
% >> proj([2, 1], [2, 2])
%
% 0 0 0 0
% 0 0 0 0
% 0 0 1 0
% 0 0 0 0
%
%
% See also BRA, KET.

rho = ket(indices, dims) * bra(indices, dims);

end
