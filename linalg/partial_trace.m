function [rho_ptr, dims_ptr] = partial_trace(rho, sys, dims)
% Return partial trace of multipartite state.
%
% Usage
% =====
%
% [RHO_PTR, DIMS_PTR] = partial_trace(RHO, SYS, DIMS)
%
% The parameter SYS contains the subsystems to trace out.
%
%
% Examples
% ========
%
% >> rho = proj([1,2,3], [3,2,3]);
% >> [rho, dims] = partial_trace(rho, [1,3], [3,2,3])
% rho =
%       0 0
%       0 1
% dims =
%       2
%
%
% >> rho = proj([1,2], [3,2]);
% >> partial_trace(rho, [1], [3,2])
% ans =
%       0 0
%       0 1
%
% >> rho = rand(16) + i*rand(16); tau = rand(16) + i*rand(16);
% >> rho = 1/trace(rho) * rho; tau = 1/trace(tau) * tau;
% >> assert_close(partial_trace(kron(rho, tau), 2, [16 16]), rho)
% >> assert_close(partial_trace(kron(rho, tau), 1, [16 16]), tau)
%
%
% See also TRACE, PARTIAL_TRANSPOSE, PARTIAL_TRACE_CHANNEL.

% reshape into bras and kets (note that the first subsystem is stored in the last component!)
rho_ptr = reshape(rho, [dims(end:-1:1), dims(end:-1:1)]);

% permute subsystems to trace over to the end
keep = [1:length(dims)];
keep(sys) = [];
keep = keep(end:-1:1);

kets_keep = length(dims) + 1 - keep;
kets_sys = length(dims) + 1 - sys;

bras_keep = length(dims) + kets_keep;
bras_sys = length(dims) + kets_sys;

pi = [kets_keep, bras_keep, kets_sys, bras_sys];
rho_ptr = permute(rho_ptr, pi);

% reshape into keep x keep x sys^2 tensor
dim_keep = prod(dims(keep));
dim_sys = prod(dims(sys));
rho_ptr = reshape(rho_ptr, [dim_keep, dim_keep, dim_sys^2]);

% contract "last two" indices
diag_indices = [1:dim_sys + 1:dim_sys^2];
rho_ptr = sum(rho_ptr(:,:,diag_indices), 3);
dims_ptr = dims(keep);

end
