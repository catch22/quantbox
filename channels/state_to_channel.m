function t_a_b = state_to_channel(rho_ab, dims_ab)
% Convert bipartite state to corresponding Choi-Jamiolkowski channel.
%
% Usage
% =====
%
% T_A_B = state_to_channel(RHO_AB, DIMS_AB)
%
%
% Examples
% ========
%
% >> assert_close(state_to_channel(maximally_entangled_state(3), [3,3]), eye(9));
%
% >> assert_close(state_to_channel(maximally_mixed_state([3, 3]), [3,3]), depolarizing_channel(3));
%
% >> d_a=3; d_b=4;
% >> for i=1:d_a; for j=1:d_b; for k=1:d_a; for l=1:d_b
% ..   t_a_b = state_to_channel(ket([i j], [d_a d_b])*bra([k l], [d_a d_b]), [d_a d_b]);
% ..   t_a_b_expected = d_a * vec(ket(j, d_b)*bra(l, d_b)) * (vec(ket(i, d_a)*bra(k, d_a))');
% ..   assert_close(t_a_b, t_a_b_expected);
% .. end; end; end; end
%
% >> d=3;
% >> for i=1:d^4
% ..   rho = zeros(d^2);
% ..   rho(i) = 1;
% ..   assert_close(channel_to_state(state_to_channel(rho, [d d])), rho);
% .. end
%
% >> d=3;
% >> for i=1:d^4
% ..   rho = zeros(d^2);
% ..   rho(i) = 1i;
% ..   assert_close(channel_to_state(state_to_channel(rho, [d,d])), rho);
% .. end
%
%
% See also CHANNEL_TO_STATE.

[dim_a, dim_b] = vunpack(dims_ab);

% illustration of the magic:
% - we start with |12><34|
% - after reshaping to [4,4,4,4], we get |2><1| at (-,-,4,3)
% - after permuting, we have |2><4| at (-,-,1,3)
% - after reshaping to [16,16], we get vec(|2><4|) at vec(|1><3|)
t_a_b = dim_a * reshape(permute(reshape(full(rho_ab), [dim_b, dim_a, dim_b, dim_a]), [1, 3, 2, 4]), [dim_b^2, dim_a^2]);

end
