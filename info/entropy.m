function S = entropy(rho)
% Return the (von Neumann) entropy of the given quantum state.
%
% Usage
% =====
%
% S = entropy(RHO)
%
%
% Examples
% ========
%
% >> entropy(proj(1, 2))
%
% ans = 0
%
% >> entropy(maximally_mixed_state(8))
%
% ans = 3
%
%
% See also MUTUAL_INFO.

% determine (non-zero) eigenvalues
p = eig(rho);
p = p(p ~= 0);

% compute entropy
S = sum(-p .* log2(p));

end
