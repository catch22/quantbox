function T = sym_to_tensor_channel(d, n)
% Return channel that extends operators on the symmetric subspace to the full tensor power.
%
% Usage
% =====
%
% T = sym_to_tensor_channel(D, N)
%
% The parameter D denotes the local dimension and N is the number of subsystems.
%
%
% Examples
% ========
%
% >> d = 2; n = 3;
% >> V = sym_to_tensor(d, n); [~, bosonic_basis_size] = vunpack(size(V));
% >> for mm = 1:bosonic_basis_size; for nn = 1:bosonic_basis_size
% ..   rho = mat(sym_to_tensor_channel(d, n) * vec(ket(mm, bosonic_basis_size) * bra(nn, bosonic_basis_size)));
% ..   rho_ = V(:,mm) * V(:,nn)';
% ..   assert_close(rho, rho_);
% .. end; end
%
%
% See also SYM_TO_TENSOR.

% get basis labels and corresponding tensors
P = partitions(d, n);
bosonic_basis_size = length(P);
bosonic_tensors = sym_to_tensor(d, n);

% fill channel column by column
tensor_basis_size = d^n;
T = zeros(tensor_basis_size^2, bosonic_basis_size ^2);
for mm = 1:bosonic_basis_size
  for nn = mm:bosonic_basis_size  % !!!
    % determine tensor product representation of |mm><nn|
    rho = bosonic_tensors(:,mm) * bosonic_tensors(:,nn)';

    % fill in columns (partial trace of |mm><nn| and of |nn><mm|)
    T(:,(nn-1)*bosonic_basis_size + mm) = vec(rho);
    T(:,(mm-1)*bosonic_basis_size + nn) = vec(rho');
  end
end

end
