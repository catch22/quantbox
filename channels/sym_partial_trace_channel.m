function T = sym_partial_trace_channel(nsys, d, n)
% Return channel corresponding to forming the partial trace for states supported on the symmetric subspace.
%
% Usage
% =====
%
% T = sym_partial_trace_channel(NSYS, D, N)
%
% The parameter NSYS denotes the number of subsystems to trace out.
% The parameter D is the local dimension and N is the total number of subsystems before tracing out.
%
%
% Examples
% ========
%
% >> T = sym_partial_trace_channel(99, 2, 100);
% >> mat(T * vec(proj(1, 101)))       % all spins up
%
% ans =
%         1 0
%         0 0
%
% >> mat(T * vec(proj(51, 101)))       % half of the spins are up
%
% ans =
%         0.50000   0.00000
%         0.00000   0.50000
%
% >> mat(T * vec(proj(101, 101)))       % all spins down
%
% ans =
%         0 0
%         0 1
%
% >> for d = 2:4; for n = 2:4; for nsys = 0:n;
% ..   dims_ab = [sym_dim(d, nsys) sym_dim(d, n - nsys)];
% ..   T = partial_trace_channel(1, dims_ab) * sym_split_channel(d, nsys, n-nsys);
% ..   assert_close(T, sym_partial_trace_channel(nsys, d, n));
% .. end; end; end
%
%
% see also PARTIAL_TRACE_CHANNEL, SYM_SPLIT_CHANNEL.

assert(nsys <= n);

% special cases
if nsys == 0
  T = speye(sym_dim(d, n)^2);
  return;
end

if nsys == n
  T = trace_channel(sym_dim(d, n));
  return
end

% fetch partitions corresponding to large and small sybspace
P_l = partitions(d, n);
P_s = partitions(d, n - nsys);
size_l = size(P_l, 1);
size_s = size(P_s, 1);

% fill channel col by col using
%
% <m| Tr(|M><N|) |n> = N * tr[ (a_1^*)^n_1 ... (a_1)^m_1 ... |M><N|]
%
% with N^(-1) = n * (n-1) * ... * (n - nsys) Prod_i Sqrt[ M_i * (M_i-1) * ... * (M_i-m_i) / m_i! ]
%             = Proc_i Sqrt[ Binom(M_i,m_i) ]
T = sparse(size_s^2, size_l^2);

for i_M = 1:size_l
  M = P_l(i_M,:);

  for i_m = 1:size_s
    mm = P_s(i_m,:);
    if any(M-mm < 0)
      continue;
    end
    for i_n = 1:size_s
      nn = P_s(i_n,:);
      N = M - mm + nn;
      [~, i_N] = ismember(N,P_l,'rows');
      assert(i_N>0);

      %val = sqrt(prod(factorial(N)) * prod(factorial(M)))/factorial(n) * factorial(nsys) / prod(factorial(M-mm)) * factorial(n-nsys)/sqrt(prod(factorial(mm))*prod(factorial(nn)));
      val = exp( sum( 0.5 * gammaln(N+1) + 0.5 * gammaln(M+1) - 0.5 * gammaln(nn+1) - 0.5 * gammaln(mm+1) - gammaln(M-mm+1) ) + gammaln(nsys+1) + gammaln(n-nsys+1) - gammaln(n+1) );
      T((i_n-1)*size_s+i_m, (i_N-1) * size_l + i_M) = val;
    end
  end
end

end
