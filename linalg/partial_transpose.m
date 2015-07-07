function rho_pt = partial_transpose(rho, sys, dims)
% Return partial transpose of multipartite state.
%
% Usage
% =====
%
% RHO_PT = partial_transpose(RHO, SYS, DIMS)
%
% The parameter SYS contains the subsystems to transpose.
%
%
% Examples
% ========
%
% >> assert_close(2 * partial_transpose(maximally_entangled_state(2), 1, [2, 2]), flip_operator([2 2]))
% >> assert_close(2 * partial_transpose(maximally_entangled_state(2), 2, [2, 2]), flip_operator([2 2]))
%
% >> rho = rand(4*3*2);
% >> assert_close(partial_transpose(rho, 1, [4, 3, 2]), transpose(partial_transpose(rho, [2 3], [4, 3, 2])))
%
%
% See also TRANSPOSE, PARTIAL_TRACE, PARTIAL_TRANSPOSE_CHANNEL.

if ~ismatrix(sys)
  sys = [sys];
end

% reshape into bras and kets (note that the first subsystem is stored in the last component!)
rho_pt = reshape(rho, [dims(end:-1:1), dims(end:-1:1)]);

% transpose by swapping the corresponding bras and kets
sys_kets = length(dims) + 1 - sys;
sys_bras = length(dims) + sys_kets;

pi = [1:2*length(dims)];
pi([sys_kets, sys_bras]) = [sys_bras, sys_kets];
rho_pt = permute(rho_pt, pi);

% reshape back
rho_pt = reshape(rho_pt, size(rho));

end
