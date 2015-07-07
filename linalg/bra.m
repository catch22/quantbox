function W = bra(indices, dims)
% Return row vector representing <indices| in computational basis.
%
% Usage
% =====
%
% W = bra(INDICES, DIMS)
%
%
% Examples
% ========
%
% >> bra([1, 2, 1], [2, 3, 2])
%
% ans = 0 0 1 0 0 0 0 0 0 0 0 0
%
% >> bra(2, 2)
%
% ans = 0 1
%
% >> bra(3, 2)
%
% ??? ...Index out of bounds...
%
%
% See also KET, PROJ.

W = ket(indices, dims)';

end
