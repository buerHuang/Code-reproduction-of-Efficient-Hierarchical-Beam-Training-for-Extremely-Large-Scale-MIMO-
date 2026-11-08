function [rate_far_hierarchy, sin_theta_best, wc_best] = ...
    training_far_hierarchy(hc, w_angle, sin_theta_grids, SNR_t, SNR_dB, Q)
% FF BTH 基线：使用每层两分支的中心子阵列角度码本。

[wc_best, sin_theta_best] = ...
    training_angle_hierarchy(hc, w_angle, sin_theta_grids, SNR_dB, Q);
rate_far_hierarchy = log2(1 + SNR_t * abs(hc * wc_best)^2);

end
