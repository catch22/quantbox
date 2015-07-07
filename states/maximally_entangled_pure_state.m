function [psi, dims] = maximally_entangled_pure_state(d)
% Return a maximally entangled state.
%
% Usage
% =====
%
% [PSI, DIMS] = maximally_entangled_pure_state(D)
%
% The parameter D denotes the local dimension.
%
%
% Examples
% ========
%
% >> [psi, dims] = maximally_entangled_pure_state(3)
%
% psi =
%       0.57735 0.00000 0.00000 0.00000 0.57735 0.00000 0.00000 0.00000 0.57735
% dims =
%       3 3
%
%
% See also MAXIMALLY_ENTANGLED_STATE.

psi = zeros(d^2, 1);
for i=1:d
  psi((i-1)*d + i) = 1;
end
psi = psi / sqrt(d);

dims = [d, d];

end
