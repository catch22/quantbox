function [rho_ext, dims_ext, rho_pt, dims_pt, info] = sym_extension(rho_ab, dims_ab, exts_ab, ppt_ab, init_ab, varargin)
% Find an extension of given bipartite density matrix to the symmetric subspace.
%
% Usage
% =====
%
% [RHO_EXT, DIMS_EXT, RHO_PT, DIMS_PT, INFO] = sym_extension(RHO_AB, DIMS_AB, EXTS_AB)
% [RHO_EXT, DIMS_EXT, RHO_PT, DIMS_PT, INFO] = sym_extension(RHO_AB, DIMS_AB, EXTS_AB, PPT_AB)
% [RHO_EXT, DIMS_EXT, RHO_PT, DIMS_PT, INFO] = sym_extension(RHO_AB, DIMS_AB, EXTS_AB, PPT_AB)
% [RHO_EXT, DIMS_EXT, RHO_PT, DIMS_PT, INFO] = sym_extension(RHO_AB, DIMS_AB, EXTS_AB, PPT_AB, INIT_AB)
% [RHO_EXT, DIMS_EXT, RHO_PT, DIMS_PT, INFO] = sym_extension(RHO_AB, DIMS_AB, EXTS_AB, PPT_AB, INIT_AB, <ARGS>)
%
% The parameter EXTS_AB specifies the number of extensions on Alice and Bobs's side.
%
% The parameter PPT_AB specifies that the partial transposition of the first PPT_AB(1) Alices
% and PPT_AB(2) Bobs should also result in a positive semidefinite matrix.
% It defaults to [0 0], which amounts to no PPT condition.
%
% The parameter INIT_AB specifies that the density operator is in fact already living in the symmetric
% subspace of INIT_AB(1) Alices of local dimension DIMS_AB(1), and likewise for Bob.
% It defaults to [1 1], which amounts to no symmetric subspace.
%
% All further arguments are passed on to SOLVE_SDP.
%
%
% Examples
% ========
%
% >> PT = channel_kron(sym_partial_trace_channel(2, 3, 3), sym_partial_trace_channel(4, 3, 5));
% >> [rho_ab, dims_ab] = choi_state(2.5);
% >> RHO_EXT = sym_extension(rho_ab, dims_ab, [3 5], [0 0], [1 1], 'sdpt3');
% ...SDPT3...
% >> assert_close(vec(rho_ab), PT*vec(RHO_EXT), 1e-08)
%
% >> RHO_EXT = sym_extension(rho_ab, dims_ab, [3 5], [0 0], [1 1], 'sedumi');
% ...SeDuMi...
% >> assert_close(vec(rho_ab), PT*vec(RHO_EXT), 1e-08)
%
% >> [rho_ab, dims_ab] = choi_state(2.5); RHO_EXT = sym_extension(rho_ab, dims_ab, [1 2], [0 1], [1 1], 'sdpt3');
% ...SDPT3...
% >> PT = channel_kron(sym_partial_trace_channel(0, 3, 1), sym_partial_trace_channel(1, 3, 2));
% >> assert_close(vec(rho_ab), PT*vec(RHO_EXT), 1e-08)
%
% >> [U, ~, ~] = svd(rand(3));
% .. sigma = U*diag([1 0 0])*U'; sigma = 1/trace(sigma)*sigma; rho_ab = kron(sigma,sigma);
% .. RHO_EXT = sym_extension(rho_ab, [3 3], [2 1], [1 0], [1 1], 'sdpt3');
% .. PT = channel_kron(sym_partial_trace_channel(1, 3, 2), sym_partial_trace_channel(0, 3, 1));
% .. assert_close(rho_ab, mat(PT*vec(RHO_EXT)), 1e-06);
% ...SDPT3...
%
% >> [U, ~, ~] = svd(rand(3) + 1i*rand(3));  % doctest: +XFAIL   % sdpt3 bug
% .. sigma = U*diag([1 0 0])*U'; sigma = 1/trace(sigma)*sigma; rho_ab = kron(sigma,sigma);
% .. RHO_EXT = sym_extension(rho_ab, [3 3], [2 1], [1 0], [1 1], 'sdpt3');
%
% >> [U, ~, ~] = svd(rand(3)); sigma = U*diag([1, 0, 0])*U'; sigma = 1/(2*trace(sigma)) * (sigma + sigma'); rho_ab_init = kron(sigma, sigma);
% >> rho_ab = sym_extension(rho_ab_init, [3 3], [2 1], [0 0], [1 1], 'sdpt3');
% ...SDPT3...
% >> RHO_EXT = sym_extension(rho_ab, [3 3], [2 4], [0 0], [2 1], 'sdpt3');
% ...SDPT3...
% >> assert_close( rho_ab_init, mat( channel_kron( sym_partial_trace_channel(1, 3, 2), sym_partial_trace_channel(3,3,4))*vec(RHO_EXT)) )
%
%
% See also SYMMETRIC_EXTENSION, SYM_PARTIAL_TRACE_CHANNEL, SYM_SPLIT_CHANNEL.

if any(init_ab > exts_ab)
  error('There must be no more initial extensions than final.')
end

% set default arguments
if nargin < 4
  ppt = [0 0];
end

if nargin < 5
  init = [1 1];
end

% determine dimensions
[dim_alice, dim_bob] = vunpack(dims_ab);
[exts_alice, exts_bob] = vunpack(exts_ab);
[ppt_alice, ppt_bob] = vunpack(ppt_ab);
[init_alice, init_bob] = vunpack(init_ab);

assert(ppt_alice <= exts_alice && ppt_bob <= exts_bob);
dims_ext = [nchoosek(dim_alice+exts_alice-1, exts_alice), nchoosek(dim_bob+exts_bob-1, exts_bob)];

% setup semidefinite program:
%
%   max 0
%   subject to
%     RHO_EXT >= 0 on Sym^exts_alice(dim_alice) ⊗ Sym^exts_bob(dim_bob)
%     tr_{2_A,...,k_A,2_B,...,k_B}(RHO_EXT) = rho_ab
%     RHO_EXT^PT >= 0   (* with respect to ppt_alice Alices and ppt_bob Bobs *)

% partial trace condition
PHI{1,1} = channel_kron(...
  sym_partial_trace_channel(exts_alice-init_alice, dim_alice, exts_alice),...
  sym_partial_trace_channel(exts_bob-init_bob, dim_bob, exts_bob));
B{1} = rho_ab;

% PPT condition
have_ppt = ppt_alice > 0 || ppt_bob > 0;
if have_ppt
  dims_pt = [nchoosek(dim_alice+ppt_alice-1, ppt_alice), nchoosek(dim_alice+exts_alice-ppt_alice-1, exts_alice-ppt_alice), nchoosek(dim_bob+ppt_bob-1, ppt_bob), nchoosek(dim_bob+exts_bob-ppt_bob-1, exts_bob-ppt_bob)];

  PHI{2,1} = -channel_kron(...
    sym_split_channel(dim_alice, ppt_alice, exts_alice - ppt_alice),...
    sym_split_channel(dim_bob, ppt_bob, exts_bob - ppt_bob));
  PHI{2,2} = channel_kron(...
    partial_transpose_channel(1, dims_pt(1:2)),...
    partial_transpose_channel(1, dims_pt(3:4)));
  B{2} = sparse(prod(dims_pt), prod(dims_pt));
end

% solve SDP
[result, info] = solve_sdp(PHI, B, [], varargin{:});

% extract result
rho_ext = result{1};
if have_ppt
  rho_pt = result{2};
end

end
