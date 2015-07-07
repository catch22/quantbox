function M = mat(V)
% Converts a vector containing the columns of a square matrix back into a matrix.
%
% Usage
% =====
%
% mat(V)
%
%
% The vector V is assumed to contain the columns of mat(V), one column after
% the other. That is, the dxd-matrix |i><j| corresponds to the vector |j,i>
% = (j-1) * d + i.
%
%
% Examples
% ========
%
% >> mat([1 2 3 4])
%
% ans =
%
%     1     3
%     2     4
%
%
% >> mat([1 2 3])
%
% ??? ...square...
%
%
% See also VEC.


dim = floor(sqrt(length(V)));
assert(dim*dim == length(V), 'Can only recover square matrices.');

M = reshape(V, dim, dim);

end
