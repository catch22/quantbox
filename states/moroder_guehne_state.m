function [erg, dims] = moroder_guehne_state(m1, m2)
% 3x3 bound entangled state (Tobias Moroder, Oleg Gittsovich, Marcus Huber and Otfried Guehne, arXiv:1405.0262).
%
% Usage
% =====
%
% [RHO, DIMS] = moroder_guehne_state(m1, m2)
%
% This state is PPT if (m1 + m2)^2 - m1*m2 <= 1, see arXiv:1405.0262.
% If the parameters do not satisfy the condition the function will fail;
% The function will assert_psd(partial_transpose( . )) before returning.
% The original numerical optimization example of the paper is obtained with:
% m1 = 0.2162; m2 = 0.4363;
%
%
% Examples
% ========
%
% This is the original numerical optimization example of the paper:
%
% >> [rho, dims] = moroder_guehne_state(0.2162,0.4363);
% >> rho_expected = [0.026, 0, 0, 0, 0.0261, 0, 0, 0, -0.0261;
% ..   0, 0.0129, 0, 0.0261, 0.0369, 0, 0, 0, 0.0369;
% ..   0, 0, 0.0129, 0, 0, -0.0369, -0.0261, 0.0369, 0;
% ..   0, 0.0261, 0, 0.0526, 0.0744, 0, 0, 0, 0.0744;
% ..   0.0261, 0.0369, 0, 0.0744, 0.132, 0, 0, 0, 0.0792;
% ..   0, 0, -0.0369, 0, 0, 0.29, 0.0744, 0.0792, 0;
% ..   0, 0, -0.0261, 0, 0, 0.0744, 0.0526, -0.0744, 0;
% ..   0, 0, 0.0369, 0, 0, 0.0792, -0.0744, 0.29, 0;
% ..  -0.0261, 0.0369, 0, 0.0744, 0.0792, 0, 0, 0, 0.132];
% >> assert_close(rho, rho_expected, 0.001);
%
%
% It is not even 2-extendible:
%
% >> [rho_ext, dims_ext, ~, ~, info] = sym_extension(rho, dims, [2 1], [0 0], [1 1], 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 1
% >> [rho_ext, dims_ext, ~, ~, info] = sym_extension(rho, dims, [2 1], [1 0], [1 1], 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 1

assert(m1>=0, 'Parameter m1 is negative.');
assert(m2>=0, 'Parameter m2 is negative.');

tmp1 = m1*m1 + m2*m2;
tmp2 = m1*m2;
m3 = sqrt( (1 - tmp1) / 2 );

assert(tmp1+tmp2<=1, 'The parameters do not satisfy the MoroderGuehne condition');

l3 = 1/ (4 + tmp2 - 2*tmp1);
l2 = 3*tmp2*l3;
l1 = 1 - l2 - 2*l3;

e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

psi1 = kron(e3,e2) + kron(e2,e3);
psi2 = kron(e1,e1) + kron(e2,e2) - kron(e3,e3);
psi3 = m1*kron(e1,e2) + m2*kron(e2,e1) + m3*(kron(e2,e2) + kron(e3,e3));
psi4 = m1*kron(e1,e3) - m2*kron(e3,e1) + m3*(kron(e3,e2) - kron(e2,e3));

erg = l1*psi1*psi1'/2 + l2*psi2*psi2'/3 + l3*psi3*psi3' + l3*psi4*psi4';
dims = [3,3];

assert_psd(partial_transpose(erg,1,dims), 1e-10, 'Ooops! Something went wrong: the matrix is not PPT when actually it is supposed to be. The code needs to be checked.');

end
