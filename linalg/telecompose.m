function [rho_ac, dims_ac] = telecompose(rho_ab, dims_ab, rho_bc, dims_bc, p_bb)
% Compose two bipartite quantum states by a teleportation-like procedure.
%
% Usage
% =====
%
% [RHO_AC, DIMS_AC] = TELECOMPOSE(RHO_AB, DIMS_AB, RHO_BC, DIMS_BC)
% [RHO_AC, DIMS_AC] = TELECOMPOSE(RHO_AB, DIMS_AB, RHO_BC, DIMS_BC, P_BB)
%
% The operator P_BB defaults to a projection onto the maximally entangled state of
% two copies of the B subsystem. It is assumed to be a projection.
%
% Return the state obtained by tensoring RHO_AB and RHO_BC, applying the projection P_BB,
% and tracing out the two B subsystems.
%
%
% Examples
% ========
%
% >> d=3;
% >> for i=1:d^4
% ..   rho = zeros(d^2); rho(i)=1;
% ..   assert_close(telecompose(maximally_entangled_state(d), [d,d], rho, [d,d]), rho/d^2);
% ..   assert_close(telecompose(rho, [d,d], maximally_entangled_state(d), [d,d]), rho/d^2);
% .. end
%
% >> d=2;
% >> for i=1:d^4; for j=1:d^4
% ..   rho_ab = zeros(d^2); rho_ab(i)=1;
% ..   rho_bc = zeros(d^2); rho_bc(j)=1;
% ..   assert_close(...
% ..     state_to_channel(telecompose(rho_ab, [d,d], rho_bc, [d,d]), [d,d]), ...
% ..     state_to_channel(rho_bc, [d,d]) * state_to_channel(rho_ab, [d,d]) / d^2);
% .. end; end
%
%
% See also TELECOMPOSE_CHANNEL.

% verify that composition is well-defined
[dim_a, dim_b1] = vunpack(dims_ab);
[dim_b2, dim_c] = vunpack(dims_bc);
assert(dim_b1 == dim_b2, 'B subsystem dimensions should agree.');

% determine projective measurement
if nargin < 5
  p_bb = maximally_entangled_state(dim_b1);
else
  assert(all(size(p_bb) == [dim_b1*dim_b2, dim_b1*dim_b2]), 'Dimension mismatch for P_BB.');
end

% project charly onto maximally entangled state (left multiplication suffices since we trace out)
rho_abbc = kron(eye(dim_a), kron(p_bb, eye(dim_c))) * kron(rho_ab, rho_bc);

% trace out Charly's subsystems
rho_ac = partial_trace(rho_abbc, [2, 3], [dims_ab, dims_bc]);
dims_ac = [dim_a, dim_c];

end
