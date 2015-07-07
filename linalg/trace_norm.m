function t = trace_norm(H)
% Return trace norm of given hermitian matrix.
%
% Usage
% =====
%
% trace_norm(H)
%
%
% The matrix H is assumed to be hermitian.
%
%
% Examples
% ========
%
% >> trace_norm([2,0;0,-9])
%
% ans = 11

t = sum(abs(eig(H)));

end
