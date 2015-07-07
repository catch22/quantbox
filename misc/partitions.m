function P = partitions(num_bins, num_items, max_items_per_bin)
% Return all partitions of given number of items into given number of bins.
%
% Usage
% =====
%
% partitions(NUM_BINS, NUM_ITEMS)
% partitions(NUM_BINS, NUM_ITEMS, MAX_ITEMS_PER_BIN)
%
%
% The parameter MAX_ITEMS_PER_BIN defaults to NUM_ITEMS.
%
%
% Examples
% ========
%
% >> partitions(2, 4)
%
% ans =
%       4 0
%       3 1
%       2 2
%       1 3
%       0 4
%
% >> partitions(4, 2)
%
% ans =
%       2 0 0 0
%       1 1 0 0
%       1 0 1 0
%       1 0 0 1
%       0 2 0 0
%       0 1 1 0
%       0 1 0 1
%       0 0 2 0
%       0 0 1 1
%       0 0 0 2
%
% >> partitions(3, 3)
%
% ans =
%       3 0 0
%       2 1 0
%       2 0 1
%       1 2 0
%       1 1 1
%       1 0 2
%       0 3 0
%       0 2 1
%       0 1 2
%       0 0 3
%
% >> partitions(3, 4, 2)
%
% ans =
%       2 2 0
%       2 1 1
%       2 0 2
%       1 2 1
%       1 1 2
%       0 2 2
%
% >> for num_items=1:3; for num_bins=1:3;
% ..   assert(length(partitions(num_bins, num_items)) == sym_dim(num_bins, num_items));
% .. end; end
%
% >> partitions(0, 4);
%
% ??? ...positive...

if nargin < 3
  max_items_per_bin = num_items;
end

assert(num_bins >= 1, 'Number of bins should be positive.');

% only one bin?
if num_bins == 1
  if num_items > max_items_per_bin
    P = [];
  else
    P = [num_items];
  end
  return
end

% recurse
P = [];
for i = min(num_items, max_items_per_bin):-1:0
  Q = partitions(num_bins-1, num_items-i, max_items_per_bin);
  [n,ans] = vunpack(size(Q));
  for j = 1:n
    P = [P; i, Q(j,:)];
  end
end

end
