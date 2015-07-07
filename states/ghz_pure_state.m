function [psi, dims] = ghz_pure_state(n, d)
% Return a multipartite GHZ state.
%
% Usage
% =====
%
% [PSI, DIMS] = ghz_pure_state(N, D)
%
% The parameter N denotes the number of subsystems.
% The parameter D denotes the local dimension.
%
%
% Examples
% ========
%
% >> [psi,dims] = ghz_pure_state(3, 2)
%
% psi =
%       0.70711 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.70711
% dims =
%       2 2 2
%
%
% See also GHZ_STATE.
dims = repmat([d], 1, n);

psi = 0;
for i=1:d
  psi = psi + ket(repmat([i], 1, n), dims);
end
psi = psi / sqrt(d);

end
