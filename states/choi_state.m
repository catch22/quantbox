function [rho, dims] = choi_state(alpha)
% Return a state from the Choi family of states.
%
% Usage
% =====
%
% [RHO, DIMS] = CHOI_STATE(ALPHA)
%
% The Choi states are a family of two-qutrit quantum states which are
%
% - separable for alpha in [2,3]
% - PPT-entangled for alpha in [1,2) and (3,4]
% - NPT for alpha in [0,1) and (4,5]
%
%
% Examples
% ========
%
% >> min(eig(partial_transpose(choi_state(0.9), 2, [3,3])))
%
% ans = -0.0029...
%
% >> min(eig(partial_transpose(choi_state(4.05), 2, [3,3])))
%
% ans = -0.0014...

% compute dimension
dims = [3 3];

% compute state
psi_plus = maximally_entangled_state(3);
sigma_plus = 1/3 * (ket([1,2], dims)*bra([1,2], dims) + ket([2,3], dims)*bra([2,3], dims) + ket([3,1], dims)*bra([3,1], dims));
V = zeros(9);
for i=1:3
  for j=1:3
    V = V + ket([i,j], dims) * bra([j,i], dims);
  end
end

rho = 2/7 * psi_plus + alpha/7 * sigma_plus + (5-alpha)/7 * V*sigma_plus*V;

end
