function [wc_best, sin_theta_best, angle_index, training_info] = ...
    training_angle_hierarchy(hc, codebook, sin_theta_grids, SNR_dB, Q)
% 按照父子索引逐层搜索角度；每一层只测量当前分支的候选波束。

layer_num = length(codebook);
angle_index = 1;
previous_num = 1;
training_info.candidate_index = cell(1, layer_num);
training_info.selected_index = zeros(1, layer_num);
training_info.overhead = 0;
training_info.feedback = layer_num;

for layer = 1:layer_num
    angle_num = size(codebook{layer}, 2);
    branch_num = angle_num / previous_num;
    candidate_index = (angle_index - 1)*branch_num + (1:branch_num);

    codebook_layer = codebook{layer}(:, candidate_index);
    power_measure = beam_measurement(hc, codebook_layer, SNR_dB, Q);
    [~, code_max] = max(power_measure);
    angle_index = candidate_index(code_max);

    training_info.candidate_index{layer} = candidate_index;
    training_info.selected_index(layer) = angle_index;
    training_info.overhead = training_info.overhead + branch_num;
    previous_num = angle_num;
end

sin_theta_best = sin_theta_grids{end}(angle_index);
wc_best = codebook{end}(:, angle_index);

end
