function result = solve_sdp_available()
% Return whether any SDP solver is available.
%
% Usage
% =====
%
% RESULT = solve_sdp_available
%
%
% Examples
% ========
%
% >> solve_sdp_available   % doctest: +XFAIL_UNLESS(solve_sdp_available)
% 1
%
% >> solve_sdp_available   % doctest: +XFAIL_IF(solve_sdp_available)
% 0
%

result = solve_sdp_sdpt3_available || solve_sdp_sedumi_available;

end
