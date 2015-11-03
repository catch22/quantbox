function [rho, dims] =  vertesi_brunner_state()
% 3x3 bound entangled state (Tamas Vertesi and Nicolas Brunner, arXiv:1405.4502).
%
% Usage
% =====
%
% [RHO, DIMS] = vertesi_brunner_state()
%
% This state is bound entangled, see arxiv:1405.4502.
%
%
% Examples
% ========
%
% This state is not even 2-extendible:
%
% >> [rho, dims] = vertesi_brunner_state();
% >> [rho_ext, dims_ext, ~, ~, info] = sym_extension(rho, dims, [2 1], [0 0], [1 1], 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 1
% >> [rho_ext, dims_ext, ~, ~, info] = sym_extension(rho, dims, [2 1], [1 0], [1 1], 'sdpt3'); info.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% ...ans = 1

a = sqrt(131/2);
l1 = 3257/6884;
l2 = 450/1721;
l3 = 450/1721;
l4 = 27/6884;

e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];

phip = kron(e1,e1) + kron(e2,e2);
phim = kron(e1,e1) - kron(e2,e2);
psip = kron(e1,e2) + kron(e2,e1);
psim = kron(e1,e2) - kron(e2,e1);
brum = 1/60*kron(e1,e3) - 3/10*kron(e3,e2);
brup = 3/10*kron(e3,e1) + 1/60*kron(e2,e3);

psi2 = a/12*psip + brum;
psi3 = a/12*phim + brup;
psi4 = kron(e3,e3) - psim;

rho = l1*phip*phip'/2 + l2*psi2*psi2' + l3*psi3*psi3' + l4*psi4*psi4'/3;
dims = [3,3];

end
