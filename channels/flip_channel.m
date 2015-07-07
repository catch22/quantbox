function T = flip_channel(dims_ab)
% Return channel that amounts to flipping the tensor factors.
%
% Usage
% =====
%
% T = flip_channels(DIMS_AB)
%
% This is the the channel T that sends |ij><kl| to |ji><lk|.
% It is equal to conjugation with the flip operator.
%
%
% Examples
% ========
%
% >> rho = rand(3*4) + 1i * rand(3*4); rho = rho + rho';
% >> rho_flipped = mat(flip_channel([3 4]) * vec(rho));
% >> assert_close(partial_trace(rho, 1, [3 4]), partial_trace(rho_flipped, 2, [4 3]));
% >> assert_close(partial_trace(rho, 2, [3 4]), partial_trace(rho_flipped, 1, [4 3]));
%
% >> rho = rand(3*4) + 1i * rand(3*4); rho = rho + rho';
% >> rho_flipped = mat(flip_channel([3 4]) * vec(rho));
% >> assert_close(rho_flipped, flip_operator([3 4]) * rho * flip_operator([3 4])');
%
%
% see also TRANSPOSE_CHANNEL

[dim_a, dim_b] = vunpack(dims_ab);

% construct appropriate permutation vector
pi = zeros(prod(dims_ab),1);
for i=1:prod(dims_ab)^2
	j4 = mod(i-1, dim_b);
	j3 = mod(floor((i-1)/dim_b), dim_a);
  j2 = mod(floor((i-1)/(dim_b*dim_a)), dim_b);
  j1=mod(floor((i-1)/(dim_b^2*dim_a)), dim_a);
	pi(i) = j2*dim_a^2*dim_b + j1*dim_a*dim_b + j4*dim_a + j3 + 1;
end

T = speye(prod(dims_ab)^2);
T = T(:,pi');

end
