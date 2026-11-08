function [power_measure, y_mean] = beam_measurement(hc, codebook, SNR_dB, Q)
% 每个波束发送 Q 个相同导频，接收端相干平均后比较功率。
% 在 fc 测量时保留原代码的 10*log10(M) 滤波增益。
% M 只在本函数内设置；噪声方差不随候选波束的接收功率改变。

assert(Q >= 1 && Q == round(Q), 'Q 必须是正整数。');
M = 1;
SNR_measure_dB = SNR_dB + 10*log10(M);
SNR = 10^(SNR_measure_dB/10);
signal = hc * codebook;
noise = sqrt(1/(2*SNR)) * ...
    (randn(Q, size(codebook, 2)) + 1j*randn(Q, size(codebook, 2)));
temp = repmat(signal, [Q, 1]) + noise;
y_mean = sum(temp, 1)/Q;
power_measure = abs(y_mean).^2;

end
