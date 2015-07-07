function [rho_ab, dims_ab] = channel_to_state(t_a_b)
% Convert channel to bipartite Choi-Jamiolkowski state.
%
% Usage
% =====
%
% [RHO_AB, DIMS_AB] = channel_to_state(T_A_B)
%
%
% The Choi-Jamiolkowski state is given by applying T_A_B to the second half of a maximally
% entangled state.
%
%
% Examples
% ========
%
% >> assert_close(channel_to_state(speye(9)), maximally_entangled_state(3));
% >> assert_close(channel_to_state(depolarizing_channel(3)), maximally_mixed_state([3 3]));
%
% >> d_a=3; d_b=4;
% >> for i=1:d_a; for j=1:d_b; for k=1:d_a; for l=1:d_b
% ..   assert_close(...
% ..     ket([i j], [d_a d_b]) * bra([k l], [d_a d_b]),...
% ..     channel_to_state(d_a * vec(ket(j, d_b) * bra(l, d_b)) * vec(ket(i,d_a) * bra(k, d_a))'));
% .. end; end; end; end
%
% >> d=3;
% >> for i=1:d^4
% ..   T = zeros(d^2);
% ..   T(i) = 1;
% ..   assert_close(state_to_channel(channel_to_state(T), [d d]), T);
% .. end
%
% >> d=3;
% >> for i=1:d^4
% ..   T = zeros(d^2);
% ..   T(i) = 1i;
% ..   assert_close(state_to_channel(channel_to_state(T), [d d]), T);
% .. end
%
%
% See also STATE_TO_CHANNEL.

% determine dimensions
[dim_b, dim_a] = vunpack(sqrt(size(t_a_b)));

% setup 1 (x) T channel
t_aa_ab = channel_kron(speye(dim_a^2), t_a_b);

% apply channel
rho_ab = mat(t_aa_ab * vec(maximally_entangled_state(dim_a)));
dims_ab = [dim_a, dim_b];

end
