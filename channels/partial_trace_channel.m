function T = partial_trace_channel(sys, dims)
% Return channel corresponding to forming the partial trace.
%
% Usage
% =====
%
% T = partial_trace_channel(SYS, DIMS)
%
%
% Examples
% ========
%
% >> full(partial_trace_channel(1, [2]))
%
% ans =
%       1 0 0 1
%
% >> rho = proj([1,2,3], [3,3,3]);
% >> partial_trace_channel([1,3], [3,3,3]) * vec(rho)
%
% ans =
%       0 0 0
%       0 1 0
%       0 0 0
%
% >> rho = proj([1,2], [3,2]);
% >> partial_trace_channel([1], [3,2]) * vec(rho)
%
% ans =
%       0 0
%       0 1
%
% >> rho=rand(16)+i*rand(16);tau=rand(16)+i*rand(16);
% >> rho = 1/trace(rho)*rho; tau = 1/trace(tau)*tau;
% >> assert_close(mat(partial_trace_channel(2,[4 4]) * vec(rho)), partial_trace(rho, 2, [4 4]))
% >> assert_close(mat(partial_trace_channel(2,[16 16]) * vec(kron(rho,tau))), rho)
% >> assert_close(partial_trace_channel(2,[4 4]) * vec(rho), channel_kron(speye(16),trace_channel(4)) * vec(rho))
%
%
% See also PARTIAL_TRACE, TRACE_CHANNEL, PARTIAL_TRANSPOSE_CHANNEL.

T = [1];
dest_dims = [];
for k=1:numel(dims)
  if ismember(k, sys)
    T = channel_kron(T, trace_channel(dims(k)));
  else
    T = channel_kron(T, eye(dims(k)^2));
    dest_dims = [dest_dims, dims(k)];
  end
end

end
