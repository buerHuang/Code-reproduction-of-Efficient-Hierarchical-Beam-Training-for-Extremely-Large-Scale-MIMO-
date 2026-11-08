clear;clc;
code_dir = fileparts(mfilename('fullpath'));
addpath(code_dir);
rng(1);

%% 系统参数：所有训练与速率计算均在 fc 处进行
Nt = 512;
fc = 100e9;
c = 3e8;
B = 0;
M = 1;
lambda_c = c/fc;
d = lambda_c/2;
Q = 16;             % 每个波束重复的导频数，沿用原 Rate_SNR 的设置
N_iter = 200;       % 正式统计时可增大到 5000
Rmin = 8;
Rmax = 30;
theta_min = -pi/3;
theta_max = pi/3;
SNR_dB_list = 15:5:35;
SNR_list = 10.^(SNR_dB_list/10);
SNR_len = length(SNR_dB_list);

%% 直接复用原文件夹的 PolarCode
s = 1;
delta = 1.2;
rho_min = 8;      % 当前 Nt、fc 下得到 16 个距离码字；距离直接由 PolarCode 决定
theta_code_min = asin(-1 + 1/Nt);
theta_code_max = asin(1 - 1/Nt);
[w_phase_near, sin_theta_near_nb, r_near_nb] = ...
    PolarCode(Nt, s, d, lambda_c, delta, rho_min, theta_code_min, theta_code_max);
w_phase_near = conj(w_phase_near);
D2 = length(sin_theta_near_nb);
S2 = length(r_near_nb);

%% EH 和 FF BTH 的中心子阵列角度码本
[w_angle_eh, sin_theta_eh, active_num_eh] = far_hierarchy_codebook(Nt, d, fc, 2);
[w_angle_bth, sin_theta_bth] = far_hierarchy_codebook(Nt, d, fc, 1);
assert(D2 == Nt && max(abs(sin_theta_eh{end} - sin_theta_near_nb)) < 1e-12, ...
    'EH 最后一层与 PolarCode 的角度网格必须一致。');
w_far = w_angle_eh{end};

%% 存储每次仿真的结果
rate_near_eh_all = zeros(SNR_len, N_iter);
rate_near_nb_all = zeros(SNR_len, N_iter);
rate_near_tp_all = zeros(SNR_len, N_iter);
rate_far_hierarchy_all = zeros(SNR_len, N_iter);
rate_opt_all = zeros(SNR_len, N_iter);
overhead_eh_all = zeros(SNR_len, N_iter);
feedback_eh_all = zeros(SNR_len, N_iter);

% 各 SNR 使用同一批用户位置，各算法使用同一个 hc。
r_list = Rmin + (Rmax - Rmin)*rand(1, N_iter);
theta_list = theta_min + (theta_max - theta_min)*rand(1, N_iter);
t0 = clock;
fprintf('Nt = %d, fc = %.1f GHz, PolarCode = %d x %d, Q = %d\n', ...
    Nt, fc/1e9, D2, S2, Q);
fprintf('EH active antennas: '); fprintf('%d ', active_num_eh); fprintf('\n');

%% SNR 循环
for idx = 1:SNR_len
    SNR_dB = SNR_dB_list(idx);
    SNR = SNR_list(idx);
    rng(1000 + idx);
    fprintf('SNR_dB = %.2f [%d/%d] | run %.2f s\n', ...
        SNR_dB, idx, SNR_len, etime(clock, t0));

    for i_iter = 1:N_iter
        r = r_list(i_iter);
        theta = theta_list(i_iter);
        [~, hc] = near_field_channel(Nt, d, fc, B, M, r, theta);

        %% EH 近场分层训练
        [rate_near_eh_all(idx, i_iter), ~, ~, ~, training_info] = ...
            training_near_eh(hc, w_angle_eh, sin_theta_eh, w_phase_near, ...
                r_near_nb, SNR, SNR_dB, Q);
        overhead_eh_all(idx, i_iter) = training_info.overhead;
        feedback_eh_all(idx, i_iter) = training_info.feedback;

        %% 近场极域码本遍历
        rate_near_nb_all(idx, i_iter) = training_near_exhaustive( ...
            hc, w_phase_near, sin_theta_near_nb, r_near_nb, SNR, SNR_dB, Q);

        %% TP：角度遍历 + 距离遍历
        rate_near_tp_all(idx, i_iter) = training_near_two_phase( ...
            hc, w_far, w_phase_near, sin_theta_near_nb, r_near_nb, SNR, SNR_dB, Q);

        %% 远场二叉分层训练
        rate_far_hierarchy_all(idx, i_iter) = training_far_hierarchy( ...
            hc, w_angle_bth, sin_theta_bth, SNR, SNR_dB, Q);

        %% 完美 CSI 下的恒模波束上界
        wc_opt = exp(-1j*angle(hc)).'/sqrt(Nt);
        rate_opt_all(idx, i_iter) = log2(1 + SNR*abs(hc*wc_opt)^2);
    end
end

%% 统计平均速率
rate_near_eh = mean(rate_near_eh_all, 2);
rate_near_nb = mean(rate_near_nb_all, 2);
rate_near_tp = mean(rate_near_tp_all, 2);
rate_far_hierarchy = mean(rate_far_hierarchy_all, 2);
rate_opt = mean(rate_opt_all, 2);
fprintf('EH training = %.2f beams, feedback = %.2f rounds, elapsed = %.2f s\n', ...
    mean(overhead_eh_all(:)), mean(feedback_eh_all(:)), etime(clock, t0));

%% 画图，沿用原 Rate_SNR 的 Plot / linspecer 写法
figure; hold on; box on; grid on;
p1 = plot(SNR_dB_list, rate_opt, 'r--', 'LineWidth', 1.6);
hold on;
p2 = plot(SNR_dB_list, rate_near_nb, 'b-s', ...
    'LineWidth', 1.6, 'MarkerFaceColor', 'w','MarkerSize',8);
p3 = plot(SNR_dB_list, rate_near_eh, 'k-o', ...
    'LineWidth', 1.6, 'MarkerFaceColor', 'w','MarkerSize',7);
p4 = plot(SNR_dB_list, rate_near_tp, 'm-^', ...
    'LineWidth', 1.6, 'MarkerFaceColor', 'w','MarkerSize',7);
p5 = plot(SNR_dB_list, rate_far_hierarchy, 'g-d', ...
    'LineWidth', 1.6, 'MarkerFaceColor', 'w','MarkerSize',7);

xlabel('SNR (dB)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Achievable Rate (bit/s/Hz)', 'Interpreter', 'latex', 'FontSize', 12);
legend([p1,p2,p3,p4,p5], ...
    {'Perfect CSI', 'NF ES', 'Proposed EH', 'NF TP', 'FF BTH'}, ...
    'Location', 'northwest', 'FontSize', 12, 'interpreter', 'latex');
set(gca, 'Position', [0.1, 0.12, 0.85, 0.85], 'FontSize', 12);
xlim([min(SNR_dB_list), max(SNR_dB_list)]);

