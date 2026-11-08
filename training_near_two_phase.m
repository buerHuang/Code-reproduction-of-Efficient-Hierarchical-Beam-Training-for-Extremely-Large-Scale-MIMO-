function [rate_near_tp, sin_theta_best, r_best, wc_best] = ...
    training_near_two_phase(hc, w_far, w_phase_near, sin_theta_grids, r_grids, SNR_t, SNR_dB, Q)
% TP 基线：全阵列远场码本遍历角度，采用3-dB中值角，再遍历该方向的全部距离码字。

%% 第一阶段：角度遍历
power_measure = beam_measurement(hc, w_far, SNR_dB, Q);

% 筛选功率不低于最大功率一半的波束索引
P_max = max(power_measure);
phi = find(power_measure >= 0.5 * P_max);

% 选取中值索引对应的角度
angle_index = round(median(phi));
sin_theta_best = sin_theta_grids(angle_index);

%% 第二阶段：距离遍历
codebook_distance = reshape(w_phase_near(:, angle_index, :), ...
    size(w_phase_near, 1), length(r_grids));
power_measure = beam_measurement(hc, codebook_distance, SNR_dB, Q);
[~, distance_index] = max(power_measure);

wc_best = codebook_distance(:, distance_index);
r_best = r_grids(distance_index) * (1 - sin_theta_best^2);
rate_near_tp = log2(1 + SNR_t * abs(hc * wc_best)^2);

end