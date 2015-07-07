function T = measurement_channel(d)
% Return channel corresponding to measurement in the computational basis.
%
% Usage
% =====
%
% T = measurement_channel(d)
%
%
% Examples
% ========
%
% >> measurement_channel(2)
%
% ans =
%        1    0   0   0
%        0    0   0   0
%        0    0   0   0
%        0    0   0   1
%
%
% See also KET, BRA, PROJ.

T = 0;
for i = 1:d
  T = T + vec(proj(i, d)) * vec(proj(i, d))';
end

end

