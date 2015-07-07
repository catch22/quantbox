function I = mutual_info(rho_ab, dims_ab)
% Return the (von Neumann) mutual information of the given bipartite quantum state.
%
% Usage
% =====
%
% I = mutual_info(RHO_AB, DIMS_AB)
%
%
% Examples
% ========
%
% >> mutual_info(maximally_entangled_state(2), [2,2])
%
% ans = 2.0000
%
% >> mutual_info(maximally_correlated_state(2), [2,2])
%
% ans = 1
%
% >> mutual_info(maximally_mixed_state([2, 2]), [2,2])
%
% ans = 0
%
%
% See also ENTROPY.

assert(length(dims_ab) == 2, 'Expect bipartite state.');

rho_a = partial_trace(rho_ab, 2, dims_ab);
rho_b = partial_trace(rho_ab, 1, dims_ab);
I = entropy(rho_a) + entropy(rho_b) - entropy(rho_ab);

end
