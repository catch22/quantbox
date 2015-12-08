function psi = random_pure_state(dims)
% Return a Haar-random pure state.
%
% Usage
% =====
%
% PSI = random_pure_state(dims)
%
%
% Examples
% ========
%
% >> psi = random_pure_state(2); size(psi), norm(psi)
% ans = 2 1
% ans = 1...
%
% >> psi = random_pure_state([2 3]); size(psi), norm(psi)
% ans = 6 1
% ans = 1...

% projector onto normalized random gaussian
D = prod(dims);
psi = normrnd(0, 1, D, 1) + 1j * normrnd(0, 1, D, 1);
psi = psi / norm(psi);

end
