function V = ket(indices, dims)
% Return column vector representing |indices> in computational basis.
%
% Usage
% =====
%
% V = ket(INDICES, DIMS)
%
%
% Examples
% ========
%
% >> ket([1, 2, 1], [2, 3, 2])
%
% ans =
%
%      0
%      0
%      1
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%      0
%
% >> ket(2, 2)
%
% ans =
%      0
%      1
%
% >> ket(3, 2)
%
% ??? ...Index out of bounds...
%
%
% See also BRA, PROJ.

assert(all(indices >= 1) && all(indices <= dims), 'Index out of bounds.');

% tensor together unit vectors
V = eye(1);
for i=1:length(indices)
  v = zeros(dims(i), 1);
  v(indices(i)) = 1;
  V = kron(V, v);
end

end
