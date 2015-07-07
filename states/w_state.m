function [rho, dims] = w_state(n)
% Return a multipartite W state.
%
% Usage
% =====
%
% [RHO, DIMS] = w_state(N)
%
% The parameter N denotes the number of subsystems.
%
%
% Examples
% ========
%
% >> [rho, dims] = w_state(3);
% >> dims
%
% dims =
%       2 2 2
%
% >> [psi, ~] = w_pure_state(3);
% >> assert_close(w_state(3), psi*psi');
%
%
% See also W_PURE_STATE.

[psi, dims] = w_pure_state(n);
rho = psi * psi';

end
