function [X, info] = solve_sdp(PHI, B, C, solver, options)
% Solve semidefinite program given in primal form.
%
% Usage
% =====
%
% [X, INFO] = solve_sdp(PHI, B)
% [X, INFO] = solve_sdp(PHI, B, C)
% [X, INFO] = solve_sdp(PHI, B, C, SOLVER)
% [X, INFO] = solve_sdp(PHI, B, C, SOLVER, OPTIONS)
%
% The semidefinite program is specified in the following primal form:
%
%   minimize    <C, X>
%   subject to  X >= 0
%               Phi(X) = B
%
% Here, <.,.> is the trace inner product.
% More generally, PHI can be a cell matrix and B, C cell arrays representing the following SDP:
%
%   minimize    sum_j <C{j}, X{j}>
%   subject to  X{j} >= 0   for all j
%               sum_{j} PHI{i,j}(X{j}) = B{i}   for all i
%
% Left-out cells in PHI default to zero superoperators.
% Left-out cells in C default to zero matrices.
% If C is not specified at all then this corresponds to a feasibility problem.
% The SOLVER parameter is either 'sdpt3' (default) or 'sedumi'.
%
% The return value INFO is solver-dependent.
%
%
% Examples
% ========
%
% >> [X, INFO] = solve_sdp(trace_channel(2), 1, diag([-1, 1]), 'sdpt3'); INFO.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% .. assert_close(X{1}, diag([1 0]), 1e-8);
% ...ans = 0
%
% >> [X, INFO] = solve_sdp(trace_channel(2), 1, diag([-1, 1]), 'sedumi');    % doctest: +SKIP_UNLESS(solve_sdp_sedumi_available)
% .. assert_close(INFO.feasratio, 1, 1e-8);
% .. assert_close(X{1}, diag([1 0]), 1e-8);
% ...
%
% >> [X, INFO] = solve_sdp(trace_channel(2), 1, complex(diag([-1, 1])), 'sdpt3'); INFO.termcode   % doctest: +SKIP_UNLESS(solve_sdp_sdpt3_available)
% .. assert_close(X{1}, diag([1 0]), 1e-8);
% ...ans = 0
%
%
% TODO
% ====
%
% - optimize constraint generation (preallocate AA and b?)
%

% determine available solvers
HAVE_SDPT3 = solve_sdp_sdpt3_available();
HAVE_SEDUMI = solve_sdp_sedumi_available();

% convert arguments to cells, if necessary, and fill in default arguments
if ~iscell(PHI);
  PHI = {PHI};
end
if ~iscell(B);
  B = {B};
end
if nargin < 5
  options = struct;
end
if nargin < 4
  if HAVE_SDPT3
    solver = 'sdpt3';
    options.printlevel = 2;
  elseif HAVE_SEDUMI
    solver = 'sedumi';
  else
    error('no solver found');
  end
end
if nargin < 3 || isempty(C)
  C = cell(1, size(PHI, 2));
elseif ~iscell(C)
  C = {C};
end

% validate sizes of cell arrays
[num_constraints, num_cones] = size(PHI);
assert(length(B) == num_constraints, 'Constraint number mismatch.');
assert(length(C) == num_cones, 'Cone number mismatch.');

% determine whether SDP is complex?
iscmp = ~(all(all(cellfun(@isreal, PHI))) && all(cellfun(@isreal, B)) && all(cellfun(@isreal, C)));

% fill in missing minimization "functionals" C{j}
for j=1:num_cones
  if isempty(C{j})
    [~, ~, sizes] = find(cellfun(@(Phi) sqrt(size(Phi, 2)), PHI(:,j)));
    assert(~isempty(sizes), sprintf('Cannot deduce size of %d-th semidefinite variable.', j));
    C{j} = sparse(sizes(1), sizes(1));
  end
  assert(size(C{j}, 1) == size(C{j}, 2), sprintf('%d-th minimization "functional" should be square matrix (got %s).', j, mat2str(size(C{j}))));
end

% fill in missing "right-hand sides" B{i}
for i=1:num_constraints
  if isempty(B{i})
    [~, ~, sizes] = find(cellfun(@(Phi) sqrt(size(Phi, 1)), PHI(i,:)));
    assert(~isempty(sizes), sprintf('Cannot deduce size of right-hand side of %d-th constraint.', i));
    B{i} = sparse(sizes(1), sizes(1));
  end
  assert(size(B{i}, 1) == size(B{i}, 2), sprintf('Right-hand side of %d-th constraint should be square matrix (got %s).', i, mat2str(size(B{i}))));
end

% fill in missing superoperator blocks PHI{i, j}
for i=1:num_constraints
  for j=1:num_cones
    if isempty(PHI{i, j})
      PHI{i, j} = sparse(numel(B{i}), numel(C{j}));
    end
  end
end

% check sizes of matrices
for i=1:num_constraints
  for j=1:num_cones
    expected = [numel(B{i}), numel(C{j})];
    got = size(PHI{i, j});
    assert(all(got == expected), sprintf('Superoperator at %s should have size %s; got %s.', mat2str([i j]), mat2str(expected), mat2str(got)));
  end
end

% construct constraints, which SDPT3 expects to be in in the form
%   sum_j <A{j,row}, X_j> = b(row)      (for all row)
% by evaluating each of the Hermitian matrix-valued constraints
%   sum_{i,j} PHI{i,j}(X{j}) = B{i}     (row all i)
% with respect to a real basis of the space of real symmetric or complex Hermitian matrices
row = 1;
for i=1:num_constraints
  t_i = size(B{i}, 1);
  for x=1:t_i
    for y=1:x
      if x == y
        % |x><x|
        idx = (x-1)*t_i + x;
        b(row) = B{i}(idx);
        for j=1:num_cones
          AA{j, row} = mat(PHI{i, j}(idx, :));
        end
      else
        % (|y><x| + |x><y|)/sqrt(2)
        idx = (y-1)*t_i + x;
        b(row) = sqrt(2) * real(B{i}(idx));
        for j=1:num_cones
          m = mat(PHI{i, j}(idx, :));
          AA{j, row} = (m + m') / sqrt(2);
        end

        % (|y><x| - |x><y|)/sqrt(2)*i
        if iscmp
          row = row + 1;
          b(row) = sqrt(2) * imag(B{i}(x, y));
          for j=1:num_cones
            m = mat(PHI{i, j}(idx, :));
            AA{j, row} = 1i * (m - m') / sqrt(2);
          end
        end
      end
      row = row + 1;
    end
  end
end
b = sparse(b);

% call solver
if strcmp(solver, 'sdpt3')
  assert(HAVE_SDPT3, 'sdpt3 not found');
  solve_sdp_with_sdpt3;
elseif strcmp(solver, 'sedumi')
  assert(HAVE_SEDUMI, 'sedumi not found');
  solve_sdp_with_sedumi;
else
  error(sprintf('unknown solver: %s', solver));
end


% ===== SDPT3 SOLVER =====
function solve_sdp_with_sdpt3

% setup semidefinite blocks
blk = cell(num_cones, 2);
for j=1:num_cones
  blk{j, 1} = 's';
  blk{j, 2} = size(C{j}, 1);
end

% sparseify right-hand side and convert to internal format expected by SDPT3
At = svec(blk, AA, ones(num_cones, 1));
clear('AA');

% call sqlp
[~, X, ~, ~, info, ~] = sqlp(blk, At, C, b, options);

% fix up trace (XXX: is this a bug in SDPT3 4.0?)
if iscmp
  for j=1:num_cones
    X{j} = 2 * X{j};
  end
end

end


% ===== SEDUMI SOLVER =====
function solve_sdp_with_sedumi

% setup semidefinite cones
assert(~iscmp, 'SeDumi solver cannot handle complex constraints (current version is buggy).');
for j=1:num_cones
  K.s(j) = size(C{j}, 1);
end

% vectorize minimization "functional"
C = cellfun(@vec, C, 'UniformOutput', false);
c = vertcat(C{:});

% build big constraint matrix
num_variables = numel(c);
num_scalar_constraints = size(AA, 2);
A = sparse(num_scalar_constraints, num_variables);
for row=1:num_scalar_constraints
  R = cellfun(@(m) vec(m), AA(:, row), 'UniformOutput', false);
  A(row, :) = vertcat(R{:})';
end
clear('AA');

% call sedumi
[x, ~, info] = sedumi(A, b, c, K, options);

% extract solution
idx = 1;
for j=1:num_cones
  X{j} = mat(x(idx:idx+numel(C{j})-1));
  idx = idx + numel(C{j});
end

end

end
