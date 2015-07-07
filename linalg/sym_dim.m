function d_sym = sym_dim(d, n)
% Return dimension of symmetric subspace.
%
% Usage
% =====
%
% D_SYM = sym_partial_trace_channel(D, N)
%
%
% Examples
% ========
%
% >> sym_dim(2, 1000)
% 1001

d_sym = nchoosek(d + n - 1, n);

end
