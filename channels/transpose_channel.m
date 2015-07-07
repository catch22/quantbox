function T = transpose_channel(d)
% Return channel corresponding to forming the transpose.
%
% Usage
% =====
%
% T = transpose_channel(d)
%
%
% Examples
% ========
%
% >> rho = rand(2);
% >> mat(transpose_channel(2)*vec(rho))-rho'
%
% ans =
%
%     0     0
%     0     0
%
% >> rho = rand(2)+i*rand(2);
% >> mat(transpose_channel(2)*vec(rho))-transpose(rho)
%
% ans =
%
%      0     0
%      0     0
%
%
% See also TRANSPOSE, PARTIAL_TRANSPOSE_CHANNEL.

% construct appropriate permutation vector
perm = zeros(d^2,1);
for i=1:d^2
	j2 = mod(i-1,d);
	j1 = mod(floor((i-1)/d),d);
	perm(i) = j2*d + j1 + 1;
end

T=speye(d^2);
T=T(perm,:);

end
