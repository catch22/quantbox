function T = partial_transpose_channel(sys, dims)
% Return channel corresponding to forming the partial transpose.
%
% Usage
% =====
%
% T = partial_transpose_channel(SYS, DIMS)
%
%
% Examples
% ========
%
% >> rho = rand(16);
% >> assert_close(mat(partial_transpose_channel(2, [4 4])*vec(rho)), partial_transpose(rho, 2, [4 4]))
%
% >> rho = rand(16) + i*rand(16);
% >> assert_close(mat(partial_transpose_channel(2, [4 4])*vec(rho)), partial_transpose(rho, 2, [4 4]))
%
%
% see also PARTIAL_TRANSPOSE, TRANSPOSE_CHANNEL, PARTIAL_TRACE_CHANNEL.

T=1;
for i=1:length(dims)
	if ismember(i, sys)
		T = channel_kron(T, transpose_channel(dims(i)));
	else
		T = channel_kron(T, speye(dims(i)^2));
	end
end
