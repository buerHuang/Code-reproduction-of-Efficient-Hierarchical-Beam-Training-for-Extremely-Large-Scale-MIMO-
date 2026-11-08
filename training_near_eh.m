function [rate_near_eh, sin_theta_best, r_best, wc_best, training_info] = ...
    training_near_eh(hc, w_angle, sin_theta_grids, w_phase_near, r_grids, SNR_t, SNR_dB, Q)
% EH：中心子阵列搜索角度，然后在原极域码本内分层搜索距离。

%% 第一阶段：角度搜索
[~, sin_theta_best, angle_index, angle_info] = ...
    training_angle_hierarchy(hc, w_angle, sin_theta_grids, SNR_dB, Q);

%% 第二阶段：对原极域码本的距离索引分组
% 码字、距离和先后次序均来自 PolarCode，不生成新的距离采样。
distance_index = 1:length(r_grids);
distance_candidates = cell(1, 0);
distance_selected = zeros(1, 0);
overhead_distance = 0;
layer = 0;

while length(distance_index) > 1
    layer = layer + 1;
    branch_num = min(4, length(distance_index));
    edge = round(linspace(0, length(distance_index), branch_num + 1));
    child_index = cell(1, branch_num);
    candidate_index = zeros(1, branch_num);

    for code = 1:branch_num
        child_index{code} = distance_index(edge(code) + 1:edge(code + 1));
        % 选取每组中间位置的现有码字作为本层代表。
        candidate_index(code) = child_index{code}(ceil(length(child_index{code})/2));
    end

    codebook_layer = reshape(w_phase_near(:, angle_index, candidate_index), ...
        size(w_phase_near, 1), branch_num);
    power_measure = beam_measurement(hc, codebook_layer, SNR_dB, Q);
    [~, code_max] = max(power_measure);
    distance_index = child_index{code_max};

    distance_candidates{layer} = candidate_index;
    distance_selected(layer) = candidate_index(code_max);
    overhead_distance = overhead_distance + branch_num;
end

%% 最终码字与速率
distance_index = distance_index(1);
wc_best = w_phase_near(:, angle_index, distance_index);
r_best = r_grids(distance_index) * (1 - sin_theta_best^2);
rate_near_eh = log2(1 + SNR_t * abs(hc * wc_best)^2);

training_info.angle_index = angle_index;
training_info.distance_index = distance_index;
training_info.angle_candidates = angle_info.candidate_index;
training_info.angle_selected = angle_info.selected_index;
training_info.distance_candidates = distance_candidates;
training_info.distance_selected = distance_selected;
training_info.overhead_angle = angle_info.overhead;
training_info.overhead_distance = overhead_distance;
training_info.overhead = angle_info.overhead + overhead_distance;
training_info.feedback = angle_info.feedback + layer;

end
