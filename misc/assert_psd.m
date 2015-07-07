function assert_psd(rho, tolerance, errmsg)
% Asserts that the given hermitian matrix is positive semidefinite up to given tolerance.
%
% Usage
% =====
%
% assert_psd(RHO)
% assert_psd(RHO, TOLERANCE)
% assert_psd(RHO, TOLERANCE, ERRMSG)
%
%
% The matrix RHO is assumed to be hermitian. The parameter TOLERANCE defaults to 1e-10.
%
%
% Examples
% ========
%
% >> assert_psd(0)
% >> assert_psd(eye(5))
%
% >> assert_psd(-eye(3))
%
% ??? ...negative...
%
% >> assert_psd(-eye(3), 1e-10, 'custom error message.')
%
% ??? ...custom...
%
%
% See also ASSERT_CLOSE.

if nargin < 2
  tolerance = 1e-10;
end

if nargin < 3
  errmsg = 'Found eigenvalues that were negative below tolerance.';
end

assert(all(eig(rho) >= -tolerance), errmsg)

end
