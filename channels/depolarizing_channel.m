function T = depolarizing_channel(d)
% Return maximally depolarizing channel in given dimension.
%
% Usage
% =====
%
% T = depolarizing_channel(D)
%
%
% Examples
% ========
%
% >> depolarizing_channel(2)
%
% ans =
%       0.50000 0.00000 0.00000 0.50000
%       0.00000 0.00000 0.00000 0.00000
%       0.00000 0.00000 0.00000 0.00000
%       0.50000 0.00000 0.00000 0.50000
%
%
% See also MAXIMALLY_MIXED_STATE.

T = vec(maximally_mixed_state(d)) * vec(eye(d))';

end

