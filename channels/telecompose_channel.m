function T = telecompose_channel(dims_ab, rho_bc, dims_bc)
% Returns channel corresponding to telecomposing with a fixed state rho_bc.
%
% Usage
% =====
%
% T = TELECOMPOSITE_CHANNEL(DIMS_AB, RHO_BC, DIMS_BC)
%
% Returns a channel T such that
%
%      mat(T * vec(RHO_AB)) == TELECOMPOSE_STATES(RHO_AB, DIMS_AB, RHO_BC, DIMS_BC)
%
% for all states RHO_AB of size DIMS_AB.
%
%
% Examples
% ========
%
% >> assert_close(telecompose_channel([3,3], maximally_entangled_state(3), [3,3]), ...
% ..   eye(81)/9)
%
% >> d=2;
% >> for i=1:d^4; for j=1:d^4
% ..   rho_ab = zeros(d^2); rho_ab(i)=1;
% ..   rho_bc = zeros(d^2); rho_bc(j)=1;
% ..   assert_close(mat(telecompose_channel([d,d], rho_bc, [d,d]) * vec(rho_ab)), ...
% ..     telecompose(rho_ab, [d,d], rho_bc, [d,d]));
% ..  end; end
%
%
% See also TELECOMPOSE.

% verify that composition is well-defined
[dim_alice, dim_b1] = vunpack(dims_ab);
[dim_b2, dim_c] = vunpack(dims_bc);
assert(dim_b1 == dim_b2);

% return channel
T = channel_kron(eye(dim_alice^2), state_to_channel(rho_bc, dims_bc) / dim_b1^2);

end
