function [rho, dims] = maximally_correlated_state(d)
% Return a maximally correlated state.
%
% Usage
% =====
%
% [RHO, DIMS] = maximally_correlated_state(D)
%
% Returns a classically maximally correlated state of the form
%
%  RHO = sum_{i=1}^d |ii><ii|.
%
%
% Examples
% ========
%
% >> maximally_correlated_state(2)
%
% ans =
%       0.50000   0.00000   0.00000   0.00000
%       0.00000   0.00000   0.00000   0.00000
%       0.00000   0.00000   0.00000   0.00000
%       0.00000   0.00000   0.00000   0.50000
%
%
% See also MAXIMALLY_ENTANGLED_STATE.

rho = 0;
for i=1:d
  rho = rho + proj([i, i], [d d]);
end
rho = rho/d;
dims = [d d];

end
