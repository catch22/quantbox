function [rho, dims] = isotropic_state(p, d)
% Return an isotropic state.
%
% Usage
% =====
%
% [RHO, DIMS] = isotropic_state(P, D)
%
% Returns a D-dimensional isotropic state according to the formula
%
%   RHO = P * maximally_entangled_state + (1 - P) * maximally_mixed_state
%
% Note: RHO is separable for a <= 1/(d + 1), and entangled otherwise.
%
%
% Examples
% ========
%
% >> assert_close(isotropic_state(1, 3), maximally_entangled_state(3));
% >> assert_close(isotropic_state(0, 3), maximally_mixed_state([3 3]));
%
%
% See also WERNER_STATE.

assert(p >= -1/(d^2-1) && p <= 1, 'Please choose p in the interval [-1/(d^2-1)..1].');

rho = p * maximally_entangled_state(d) + (1-p) * maximally_mixed_state([d, d]);
dims = [d, d];

end
