function key = devetak_winter(rho_abe, dims_abe, lambda_a_x)
% Return the Devetak-Winter lower bound on the one-way secret key rate.
%
% Usage
% =====
%
% KEY = devetak_winter(RHO_ABE, DIMS_ABE, LAMBDA_A_X)
%
% The parameter LAMBDA_A_X is assumed to be a quantum-to-classical map
% from A to X. The Devetak-Winter bound is computed for the cqq state
% resulting from applying this map to RHO_ABE.
%
% See Devetak and Winter (2003) for details.
%
% Examples
% ========
%
% >> devetak_winter(maximally_mixed_state([2,2,2]), [2,2,2], speye(2^2))
%
% ans = 0
%
% >> devetak_winter(kron(maximally_correlated_state(2), maximally_mixed_state(2)), [2,2,2], speye(2^2))
%
% ans = 1
%
% >> devetak_winter(kron(maximally_entangled_state(2), proj(1, 2)), [2,2,2], measurement_channel(2))
%
% ans = 1.00000
%
% >> devetak_winter(kron(maximally_entangled_state(4), proj(1, 3)), [4,4,3], measurement_channel(4))
%
% ans = 2

% determine dimensions
[dim_a, dim_b, dim_e] = vunpack(dims_abe);
assert_close(dim_a*dim_a, size(lambda_a_x, 2));
dim_x = sqrt(size(lambda_a_x, 1));

% compute cqq state
lambda_abe_xbe = channel_kron(lambda_a_x, speye(dim_b * dim_b * dim_e * dim_e));
dims_xbe = [dim_x, dim_b, dim_e];
rho_xbe = mat(lambda_abe_xbe * vec(rho_abe));

% compute partial traces
dims_xb = dims_xbe([1 2]);
rho_xb = partial_trace(rho_xbe, 3, dims_xbe);

dims_xe = dims_xbe([1 3]);
rho_xe = partial_trace(rho_xbe, 2, dims_xbe);

% compute Devetak-Winter lower bound
key = mutual_info(rho_xb, dims_xb) - mutual_info(rho_xe, dims_xe);

end
