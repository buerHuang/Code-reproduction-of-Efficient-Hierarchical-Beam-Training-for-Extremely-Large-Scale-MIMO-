function [codebook, sin_theta_grids, active_num] = far_hierarchy_codebook(Nt, d, fc, layer_step)
% 中心子阵列角度码本。layer_step = 2 为 EH，layer_step = 1 为二叉分层。

P = round(log2(Nt));
assert(Nt >= 2 && 2^P == Nt, 'Nt 必须是大于等于 2 的二次幂。');
assert(layer_step == 1 || layer_step == 2, 'layer_step 只能取 1 或 2。');

layer_num = ceil(P / layer_step);
codebook = cell(1, layer_num);
sin_theta_grids = cell(1, layer_num);
active_num = zeros(1, layer_num);
nn = (-(Nt - 1)/2:(Nt - 1)/2)';
c = 3e8;

for layer = 1:layer_num
    active_num(layer) = 2^min(layer_step * layer, P);
    K = active_num(layer);
    sin_theta_grids{layer} = -1 + (2*(1:K) - 1)/K;

    active_index = (Nt - K)/2 + (1:K);
    codebook{layer} = zeros(Nt, K);
    % hc 是行信道，发射码字使用与信道相反的相位。
    codebook{layer}(active_index, :) = ...
        exp(-1j*2*pi*fc/c*d*nn(active_index)*sin_theta_grids{layer})/sqrt(K);
end

end
