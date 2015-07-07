function T = trace_channel(d)
% Return channel corresponding to forming the trace.
%
% Usage
% =====
%
% T = trace_channel(d)
%
%
% Examples
% ========
%
% >> trace_channel(2)
%
% ans =
%       1 0 0 1
%
% >> trace_channel(3) * vec([1 2 3; 4 5 6; 7 8 9])
%
% ans = 15
%
%
% See also TRACE, PARTIAL_TRACE_CHANNEL.

T = vec(eye(d))';

end
