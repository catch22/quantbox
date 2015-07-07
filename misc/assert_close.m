function assert_close(a, b, tolerance, errmsg)
% Asserts that the given parameters are equal up to given tolerance.
%
% Usage
% =====
%
% assert_close(A, B)
% assert_close(A, B, TOLERANCE)
% assert_close(A, B, TOLERANCE, ERRMSG)
%
%
% The parameter TOLERANCE defaults to 1e-10.
%
%
% Examples
% ========
%
% >> assert_close(0, 1e-5, 1e-4)
%
% >> assert_close(0, 1)
%
% ??? ...differ...
%
% >> assert_close(0, 1, 1e-10, 'custom error message.')
%
% ??? ...custom...
%
%
% See also ASSERT_PSD.

if nargin < 3
  tolerance = 1e-10;
end

if nargin < 4
  errmsg = 'Inputs differ more than tolerance.';
end

delta = reshape(abs(a - b), [], 1);
assert(normest(delta, tolerance / 2) <= tolerance, errmsg);

end
