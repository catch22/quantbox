function [rho, dims] = werner_state(p, d)
% Return a Werner state.
%
% Usage
% =====
%
% [RHO, DIMS] = werner_state(P, D)
%
% Returns a D-dimensional Werner state according to the formula
%
%     RHO = p / (D^2 + D) * (id + F) + (1 - p) / (D^2 - D) * (id - F)
%
% where F denotes the flip operator.
%
% Examples
% ========
%
% >> werner_state(0, 2)
%
% ans =
%         0.00000  0.00000  0.00000  0.00000
%         0.00000  0.50000 -0.50000  0.00000
%         0.00000 -0.50000  0.50000  0.00000
%         0.00000  0.00000  0.00000  0.00000
%
%
% See also ISOTROPIC_STATE.


% assemble state
F = flip_operator([d, d]);
rho = p * (1/(d^2 + d)) * (eye(d^2) + F) + (1-p) * (1/(d^2 - d)) * (eye(d * d) - F);
dims = [d, d];

end
