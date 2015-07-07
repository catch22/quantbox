function F = flip_operator(dims_ab)
% Return the flip operator.
%
% Usage
% =====
%
% F = flip_operator(DIMS_AB)
%
% Returns the unitary F that sends |ij> to |ji>.
%
%
% Examples
% ========
%
% >> flip_operator([2 2])
% ans =
%        1   0   0   0
%        0   0   1   0
%        0   1   0   0
%        0   0   0   1

[dim_a, dim_b] = vunpack(dims_ab);
dims_ba = [dim_b dim_a];

F = 0;
for i=1:dim_a
  for j=1:dim_b
    F = F + ket([j i], dims_ba) * bra([i j], dims_ab);
  end
end

end
