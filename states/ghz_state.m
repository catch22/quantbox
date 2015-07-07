function [rho, dims] = ghz_state(n, d)
% Return a multipartite GHZ state.
%
% Usage
% =====
%
% [RHO, DIMS] = ghz_state(N, D)
%
% The parameter N denotes the number of subsystems.
% The parameter D denotes the local dimension.
%
%
% Examples
% ========
%
% >> [rho, dims] = ghz_state(3, 2)
%
% rho =
%       0.50000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.50000
%       0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
%       0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
%       0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
%       0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
%       0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
%       0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
%       0.50000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.50000
% dims =
%       2 2 2
%
%
% See also GHZ_PURE_STATE.

[psi, dims] = ghz_pure_state(n, d);
rho = psi * psi';

end
