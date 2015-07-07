function varargout = vunpack(what)
% Unpack matrix entries or cell array elements into variables.
%
% Usage
% =====
%
% [X, Y, ..] = vunpack(WHAT)
%
%
% The number of return values has to be equal to the number of elements in WHAT.
%
%
% Examples
% ========
%
% >> [a,b,c] = vunpack({1,3,5})
%
% a = 1
% b = 3
% c = 5
%
% >> [a,b,c] = vunpack([2,4,6])
%
% a = 2
% b = 4
% c = 6
%
% >> a = vunpack(12)
%
% a = 12
%
% >> vunpack([2, 4])
%
% ??? ...Number of return values...

assert(nargout == length(what), 'Number of return values does not match number of elements in input.');

if iscell(what(1))
  for i=1:length(what)
    varargout{i} = what{i};
  end
else
  for i=1:length(what)
    varargout{i} = what(i);
  end
end

end
