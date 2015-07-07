function [rho, dims] = maximally_entangled_state(d)
% Return a maximally entangled state.
%
% Usage
% =====
%
% [RHO, DIMS] = maximally_entangled_state(D)
%
% The parameter D denotes the local dimension.
%
%
% Examples
% ========
%
% >> [rho, dims] = maximally_entangled_state(2)
%
% rho =
%       0.50000 0.00000 0.00000 0.50000
%       0.00000 0.00000 0.00000 0.00000
%       0.00000 0.00000 0.00000 0.00000
%       0.50000 0.00000 0.00000 0.50000
% dims =
%       2 2
%
%
% See also MAXIMALLY_ENTANGLED_PURE_STATE, MAXIMALLY_CORRELATED_STATE.


[psi, dims] = maximally_entangled_pure_state(d);
rho = psi * psi';

end
