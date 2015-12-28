function I = coherent_info(rho_ab, dims_ab)
%% Return the coherent information I(A>B) of the given bipartite quantum state.
%
% Usage
% =====
%
% I = coherent_info(RHO_AB, DIMS_AB)
%
%
% Examples
% ========
%
% >> coherent_info(maximally_entangled_state(2), [2,2])
%
% ans = 1...
%
% >> coherent_info(maximally_correlated_state(2), [2,2])
%
% ans = 0
%
% >> coherent_info(maximally_mixed_state([2, 2]), [2,2])
%
% ans = -1
%
%
% See also ENTROPY and MUTUAL_INFO.

assert(length(dims_ab) == 2, 'Expect bipartite state.');

rho_b = partial_trace(rho_ab, 1, dims_ab);
I = entropy(rho_b) - entropy(rho_ab);

end
