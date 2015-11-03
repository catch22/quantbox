function [rho_ext, dims_ext, info] = permsym_extension(rho_ab, dims_ab, k, varargin)
% Find permutation-symmetric k-extension of given bipartite density matrix.
%
% Usage
% =====
%
% [RHO_EXT, DIMS_EXT, INFO] = permsym_extension(RHO_AB, DIMS_AB, K)
% [RHO_EXT, DIMS_EXT, INFO] = permsym_extension(RHO_AB, DIMS_AB, K, <ARGS>)
%
% All further arguments are passed on directly to SOLVE_SDP.
%
%
% Examples
% ========
%
% >> [rho_ext, ~, info] = permsym_extension(proj([1,2], [2,2]), [2,2], 2, 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 0
% >> assert_close(rho_ext, proj([1,2,2], [2,2,2]), 1e-8);                                           % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
%
% >> [rho_ext, ~, info] = permsym_extension(proj([1,2], [2,2]), [2,2], 2, 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 0
% >> assert_close(rho_ext, proj([1,2,2], [2,2,2]), 1e-8);                                           % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
%
% >> [~, ~, info] = permsym_extension(maximally_entangled_state(2), [2,2], 2, 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 1
%
%
% See also BOSONIC_EXTENSION.

% determine source and target dimensions
[dim_a, dim_b] = vunpack(dims_ab);
dims_ext = [dim_a, repmat(dim_b, [1 k])];

% set up the following semidefinite program:
%
%   max
%     0
%   subject to
%     partially_trace_out_all_bobs_except_j(rho_ext) = rho  (for all j=1..k)
%     rho_ext >= 0

for i = 1:k
  PHI{i,1} = partial_trace_channel(1 + [1:(i-1), (i+1):k], dims_ext);
  B{i} = rho_ab;
end

% run sedumi and extract result
[result, info] = solve_sdp(PHI, B, [], varargin{:});
rho_ext = result{1};

% symmetrize with respect to the Bobs
X = zeros(prod(dims_ext));
PI = perms(1:k);
for i=1:length(PI)
  pi = PI(i, :);
  Y = vec(rho_ext);
  Y = reshape(Y, [dims_ext dims_ext]);
  Y = permute(Y, [pi, k+1, k+1+pi, 2*(k+1)]);
  Y = reshape(Y, [prod(dims_ext) prod(dims_ext)]);
  X = X + Y;
end
rho_ext = X / length(PI);

end
