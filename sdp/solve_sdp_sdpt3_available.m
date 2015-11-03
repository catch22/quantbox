function result = solve_sdp_sdpt3_available()
% Return whether the SDPT3 SDP solver is available.
%
% Usage
% =====
%
% RESULT = solve_sdp_sdpt3_available
%
%
% Examples
% ========
%
% >> solve_sdp_sdpt3_available   % doctest: +XFAIL_UNLESS(solve_sdp_sdpt3_available)
% 1
%
% >> solve_sdp_sdpt3_available   % doctest: +XFAIL_IF(solve_sdp_sdpt3_available)
% 0
%

result = exist('sdpt3', 'file') ~= 0;

end
