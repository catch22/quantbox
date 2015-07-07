function T_ab_cd = channel_kron(T_a_c, T_b_d)
% Return tensor product of two channels.
%
% Usage
% =====
%
% T_AB_CD = channel_kron(T_A_C, T_B_D)
%
% Note that KRON cannot be used directly as it is not compatible with the bases for vectorized density matrices.
%
%
% Examples
% ========
%
% >> rho = rand(16)+i*rand(16);tau=rand(16)+i*rand(16);
% >> assert_close(mat(channel_kron(speye(16), transpose_channel(4)) * vec(rho)), partial_transpose(rho, 2, [4 4]))
% >> assert_close(mat(channel_kron(speye(256), transpose_channel(16)) * vec(kron(rho,tau))), kron(rho, transpose(tau)))
% >> assert_close(mat(channel_kron(speye(16), trace_channel(4)) * vec(rho)), partial_trace(rho, 2, [4 4]))
% >> assert_close(mat(channel_kron(speye(256), trace_channel(16)) * vec(kron(rho,tau))), rho * trace(tau))

% determine dimensions of source and range of channels
[dims_c, dims_a] = vunpack(sqrt(size(T_a_c)));
[dims_d, dims_b] = vunpack(sqrt(size(T_b_d)));

% take kronecker product and change back to correct basis
T_ab_cd = basis_change([dims_c, dims_d]) * kron(T_a_c, T_b_d) * basis_change([dims_a dims_b])';

end


function U=basis_change(dimsOUT)
  perm_vec=zeros(prod(dimsOUT)^2,1);

  for i=1:prod(dimsOUT)^2
    %% decompose index i into corresponding mixed number system
    % conversion form kron basis to comp basis
    j4=mod(i-1,dimsOUT(2));
    j3=mod(floor((i-1)/dimsOUT(2)),dimsOUT(2));
    j2=mod(floor((i-1)/(dimsOUT(2)^2)),dimsOUT(1));
    j1=mod(floor((i-1)/(dimsOUT(2)^2*dimsOUT(1))),dimsOUT(1));
    perm_vec(i) = j1*(dimsOUT(2)^2*dimsOUT(1))+j3*prod(dimsOUT)+j2*dimsOUT(2)+j4+1;
    %conv_op_OUT(j1*(dimsOUT(2)^2*dimsOUT(1))+j3*prod(dimsOUT)+j2*dimsOUT(2)+j4+1,i)=1;
  end
  U = speye(prod(dimsOUT)^2);
  U = U(:,perm_vec');
end
