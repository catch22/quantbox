function T = sym_split_channel(d, k, l)
% Embed symmetric subspace into tensor product of symmetric subspaces.
%
% Usage
% =====
%
% T = sym_split_channel(D, K, L)
%
% The parameter D is the local dimension.
%
% Returns channel that isometrically embeds Sym^(K+L)(C^D) into
% Sym^K(C^D) tensor Sym^L(C^D).
%
%
% Examples
% ========
%
% >> T = sym_split_channel(2, 1, 2);
% >> assert(isequal(sqrt(size(T)), [sym_dim(2, 1) * sym_dim(2, 2), sym_dim(2, 3)]));
% >> assert_close(sym_split_channel(3,1,1),sym_to_tensor_channel(3,2))
%
% >> [U,~,~] = svd(rand(2)+i*rand(2)); % test random state with complex entries
% >> sigma = U*diag([1,0])*U'; rho = kron(sigma,sigma); % generate bipartite product state
% >> RHO = mat(channel_kron(speye(4), sym_split_channel(2,1,1)') * vec(kron(rho,sigma))); % ... and a possible extension on Bob's part
% >> assert_close(mat(channel_kron(speye(4), sym_partial_trace_channel(1,2,2))*vec(RHO)),rho, 1e-9)
%
% >> rho = sym_split_channel(3, 5, 10);  % speed test
%
% >> rho = sym_split_channel(2, 1, 199);
% >> assert(~any(any(isinf(rho))));
% >> assert(~any(any(isnan(rho))));
%
%
% See also SYM_PARTIAL_TRACE_CHANNEL.

% generate occupation number basis labels for all systems
P_k = partitions(d, k);
P_l = partitions(d, l);
P = partitions(d, k + l);
basis_size_k = length(P_k(:, 1));
basis_size_l = length(P_l(:, 1));
basis_size = length(P(:, 1));

% generate occupation number basis labels for all systems
P_k = partitions(d, k);
P_l = partitions(d, l);
P = partitions(d, k + l);
basis_size_k = length(P_k(:, 1));
basis_size_l = length(P_l(:, 1));
basis_size = length(P(:, 1));

% nothing to do?
if k==0 || l==0
  T = speye(basis_size^2);
  return
end

% build basis of occupation numbers in tensor product space
T = sparse(basis_size_k*basis_size_l, basis_size);
for i_k=1:basis_size_k
  p_k = P_k(i_k, :);
  for i_l=1:basis_size_l
    p_l = P_l(i_l, :);
    p = p_k + p_l;
    [~, i] = ismember(p, P, 'rows');

    %val = sqrt(factorial(k)*factorial(l) / (prod(factorial(p_k))*prod(factorial(p_l))) * prod(factorial(p)) / factorial(k+l));
    val = exp(0.5*(gammaln(k+1) + gammaln(l+1) - sum(gammaln(p_k+1)) - sum(gammaln(p_l+1)) ...
      + sum(gammaln(p+1)) - gammaln(k+l+1)));
    T(basis_size_l * (i_k - 1) + i_l, i) = val;
  end
end


T = kron(T, T);

end
