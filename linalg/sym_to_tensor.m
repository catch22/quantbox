function V = sym_to_tensor(d, n)
% Return isometric embedding of symmetric subspace into the full tensor power.
%
% Usage
% =====
%
% V = sym_to_tensor(D, N)
%
% The parameter D denotes the local dimension and N is the number of subsystems.
%
% The columns of V contain the coefficients of the canonical occupation number states
% in the tensor product basis.
%
% Examples
% ========
%
% >> for d = [2,3]; for n = [1,3]
% ..   V = sym_to_tensor(d, n);
% ..   assert_close(V'*V, eye(sym_dim(d, n)));  % should better be an isometry
% .. end; end
%
% >> V = sym_to_tensor(2, 3); dims = [2,2,2];
% >> assert_close(V(:, 1), ket([1,1,1], dims));   % |3,0>
% >> assert_close(V(:, 2), 1/sqrt(3)*(ket([1,1,2], dims) + ket([1,2,1], dims) + ket([2,1,1], dims)));  % |2,0>
%
%
% See also SYM_TO_TENSOR_CHANNEL.

% get bosonic basis
P = partitions(d, n);
P_size = length(P);

% compute corresponding tensors, i.e. tensor<-bosonic basis change matrix
tensor_basis_size = d^n;
V = zeros(tensor_basis_size, P_size);
PI = perms(1:n);
dims = repmat(d, [1,n]);
for m = 1:P_size
  M = P(m,:);

  % set indices to (1^{m_1} ... d^{m_d})
  indices = [];
  for i = 1:d
    indices = [indices, repmat(i, [1,M(i)])];
  end

  % sum up all permutations
  v = zeros(tensor_basis_size, 1);
  for i = 1:length(PI)
    v = v + ket(indices(PI(i,:)), dims);
  end

  % normalize and store in matrix
  v = v / sqrt(prod(factorial(M)) * factorial(n));
  V(:,m) = v;
end

end
