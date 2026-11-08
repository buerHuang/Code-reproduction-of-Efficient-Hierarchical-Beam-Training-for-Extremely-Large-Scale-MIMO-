function [rate_near_nb, sin_theta_best, r_best, wc_best] = ...
    training_near_exhaustive(hc, w_phase_near, sin_theta_grids, r_grids, SNR_t, SNR_dB, Q)
% 在同一份 PolarCode 极域码本中遍历全部角度和距离。

Nt = size(w_phase_near, 1);
codebook = reshape(w_phase_near, Nt, []);
power_measure = beam_measurement(hc, codebook, SNR_dB, Q);
[~, code_max] = max(power_measure);
[angle_index, distance_index] = ind2sub( ...
    [length(sin_theta_grids), length(r_grids)], code_max);

wc_best = codebook(:, code_max);
sin_theta_best = sin_theta_grids(angle_index);
r_best = r_grids(distance_index) * (1 - sin_theta_best^2);
rate_near_nb = log2(1 + SNR_t * abs(hc * wc_best)^2);

end
