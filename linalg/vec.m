function V = vec(M)
% Converts a square matrix to a vector containing the matrix' columns.
%
% Usage
% =====
%
% vec(M)
%
%
% The matrix M is assumed to be a square matrix. The resulting vector stores
% the matrix elements column by column, in agreement with Octave/Matlab's
% conventions. That is, the dxd-matrix |i><j| corresponds to the vector |j,i>
% = (j-1) * d + i.
%
%
% Examples
% ========
%
% >> rho = [[1;2], [3;4]]; vec(rho)
% ans =
%
%      1
%      2
%      3
%      4
%
%
% See also MAT.

V = reshape(M, [], 1);

end
