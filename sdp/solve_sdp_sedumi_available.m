function result = solve_sdp_sedumi_available()
% Return whether the SeDumi SDP solver is available.
%
% Usage
% =====
%
% RESULT = solve_sdp_sedumi_available
%
%
% Examples
% ========
%
% >> solve_sdp_sedumi_available   % doctest: +XFAIL_UNLESS(solve_sdp_sedumi_available)
% 1
%
% >> solve_sdp_sedumi_available   % doctest: +XFAIL_IF(solve_sdp_sedumi_available)
% 0
%

result = exist('sedumi', 'file') ~= 0;

end
