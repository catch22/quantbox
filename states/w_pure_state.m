function [psi, dims] = w_pure_state(n)
% Return a multipartite W state.
%
% Usage
% =====
%
% [RHO, DIMS] = w_pure_state(N)
%
% The parameter N denotes the number of subsystems.
%
%
% Examples
% ========
%
% >> [psi,dims] = w_pure_state(3)
%
% psi =
%       0.00000 0.57735 0.57735 0.00000 0.57735 0.00000 0.00000 0.00000
% dims =
%       2 2 2
%
% >> [psi,dims] = w_pure_state(4)
%
% psi =
%       0.00000 0.50000 0.50000 0.00000 0.50000 0.00000 0.00000 0.00000 0.50000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000
% dims =
%       2 2 2 2
%
%
% See also W_STATE.

dims = repmat([2], 1, n);

psi = 0;
for i=1:n
  indices = ones(1, n);
  indices(i) = 2;
  psi = psi + ket(indices, dims);
end
psi = psi / sqrt(n);

end
