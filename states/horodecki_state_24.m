function [rho, dims] = horodecki_state_24(b)
% 2x4 bound entangled state (P. Horodecki, 1997).
%
% Usage
% =====
%
% [RHO, DIMS] = horodecki_state_24(b)
%
% This state is bound entangled for 0 < a < 1, cf. (1.27) in [Bäuml].
%
%
% Examples
% ========
%
% >> [rho, dims] = horodecki_state_24(0.5);
% >> [rho_ext, dims_ext, ~, ~, info] = sym_extension(rho, dims, [11 1], [0 0], [1 1], 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 0
% >> [rho_ext, dims_ext, ~, ~, info] = sym_extension(rho, dims, [12 1], [0 0], [1 1], 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 1
% >> [rho_ext, dims_ext, ~, ~, info] = sym_extension(rho, dims, [2 1], [1 0], [1 1], 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 1

e1 = [1;0];
e2 = [0;1];
q1 = [1;0;0;0];
q2 = [0;1;0;0];
q3 = [0;0;1;0];
q4 = [0;0;0;1];

Psi1 = kron(e1, q1) + kron(e2, q2);
Psi2 = kron(e1, q2) + kron(e2, q3);
Psi3 = kron(e1, q3) + kron(e2, q4);
q14 = kron(e1, q4);

Phi_b = kron(e2, sqrt((1 + b) / 2) * q1 + sqrt((1-b) / 2) * q4);

sigma_insep = (Psi1 * Psi1' + Psi2 * Psi2' + Psi3 * Psi3' + q14 * q14');

rho = b * sigma_insep / (7 * b + 1) + Phi_b * Phi_b' / (7 * b + 1);

dims = [2,4];

end
