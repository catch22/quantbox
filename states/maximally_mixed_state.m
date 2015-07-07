function [rho, dims] = maximally_mixed_state(dims)
% Return a maximally mixed state.
%
% Usage
% =====
%
% [RHO, DIMS] = maximally_mixed_state(DIMS)
%
%
% Examples
% ========
%
% >> [rho, dims] = maximally_mixed_state([2, 2]); full(rho), dims
%
% ans =
%       0.25000 0.00000 0.00000 0.00000
%       0.00000 0.25000 0.00000 0.00000
%       0.00000 0.00000 0.25000 0.00000
%       0.00000 0.00000 0.00000 0.25000
% dims =
%       2 2

d = prod(dims);
rho = eye(d) / d;

end
