%% =====================================================
%  主程序 - LFM无线电引信仿真
%  =====================================================
% 完整的仿真流程：
% 1. 加载参数
% 2. 生成LFM信号和目标回波
% 3. 加入干扰和噪声
% 4. 匹配滤波和脉冲压缩
% 5. 距离-多普勒处理
% 6. 引信判决
% 7. 绘制结果

clear all; close all; clc;

%% ===== 第1步：加载参数 =====
parameter;

%% ===== 第2步：生成LFM信号 =====
fprintf('生成LFM信号...\n');
t = (0:Ns-1) * dt;
s_tx = generateLFM(f0, B, T, Tr, Np, fs, Ns);

if plot_LFM
    figure('Name','LFM信号');
    subplot(2,2,1);
    t_plot = t(1:min(length(t), 200))*1e6;
    s_plot = s_tx(1:min(length(s_tx), 200));
    plot(t_plot, real(s_plot));
    xlabel('时间 (μs)'); ylabel('幅度');
    title('LFM信号 - 实部');
    grid on;
    
    subplot(2,2,2);
    plot(t_plot, imag(s_plot));
    xlabel('时间 (μs)'); ylabel('幅度');
    title('LFM信号 - 虚部');
    grid on;
    
    subplot(2,2,3);
    plot(t_plot, abs(s_plot));
    xlabel('时间 (μs)'); ylabel('幅度');
    title('LFM信号 - 包络');
    grid on;
    
    subplot(2,2,4);
    [~, f_plot, ~, P_plot] = spectrogram(s_plot, 32, 16, 128, fs);
    imagesc(t_plot, f_plot/1e9, 10*log10(P_plot + eps));
    colorbar;
    xlabel('时间 (μs)'); ylabel('频率 (GHz)');
    title('LFM信号 - 频谱图');
    axis xy;
end

%% ===== 第3步：生成目标回波 =====
fprintf('生成目标回波...\n');
r_echo = targetEcho(s_tx, R_target, v_target, f0, c, sigma_rcs);

if plot_Echo
    figure('Name','目标回波');
    subplot(1,2,1);
    t_zoom = min(3e-6, t(end));  % 放大前3μs
    idx_zoom = find(t <= t_zoom);
    plot(t(idx_zoom)*1e6, real(r_echo(idx_zoom)));
    xlabel('时间 (μs)'); ylabel('幅度');
    title('目标回波 - 实部');
    grid on;
    
    subplot(1,2,2);
    plot(t(idx_zoom)*1e6, abs(r_echo(idx_zoom)));
    xlabel('时间 (μs)'); ylabel('幅度');
    title('目标回波 - 包络');
    grid on;
end

%% ===== 第4步：生成干扰信号 =====
fprintf('生成干扰信号...\n');
i_total = zeros(1, Ns);

if enable_DRFM
    fprintf('  DRFM干扰...\n');
    i_drfm = DRFM(r_echo, drfm_delay, JNR_DRFM, fs, drfm_phase_offset);
    i_total = i_total + i_drfm;
end

if enable_RGPO
    fprintf('  RGPO干扰...\n');
    i_rgpo = RGPO(R_target, n_gates_RGPO, gate_spacing, JNR_RGPO, f0, c, fs, Np, Tr);
    i_total = i_total + i_rgpo;
end

if enable_VGPO
    fprintf('  VGPO干扰...\n');
    i_vgpo = VGPO(s_tx, R_target, n_doppler_VGPO, JNR_VGPO, f0, c, Tr, Np);
    i_total = i_total + i_vgpo;
end

%% ===== 第5步：加入噪声 =====
fprintf('加入高斯白噪声...\n');
n_noise = addNoise(Ns, sigma_n);

%% ===== 第6步：合成接收信号 =====
fprintf('合成接收信号...\n');
y_rx = r_echo + i_total + n_noise;

%% ===== 第7步：匹配滤波和脉冲压缩 =====
fprintf('匹配滤波和脉冲压缩...\n');
y_compressed = matchedFilter(y_rx, s_tx, fs);

if plot_RangeProfile
    figure('Name','距离轮廓');
    range_plot = (0:length(y_compressed)-1) * c / (2*fs);
    plot(range_plot, 20*log10(abs(y_compressed) + eps));
    xlabel('距离 (m)'); ylabel('幅度 (dB)');
    title('距离轮廓（压缩后）');
    grid on;
    xlim([0 1000]);
end

%% ===== 第8步：距离-多普勒处理 =====
fprintf('距离-多普勒处理...\n');
rd_map = rangeDoppler(y_compressed, Np, Tr, Nfft_range, Nfft_doppler);

if plot_RDMap
    figure('Name','距离-多普勒图');
    range_axis_plot = (0:Nfft_range-1) * c / (2*fs);
    velocity_axis_plot = (0:Nfft_doppler-1) / (Nfft_doppler * Tr) * (c / (2*f0));
    
    imagesc(velocity_axis_plot, range_axis_plot, 20*log10(rd_map + eps));
    colorbar;
    xlabel('速度 (m/s)'); ylabel('距离 (m)');
    title('距离-多普勒图');
    axis xy;
    set(gca, 'YLim', [0 1000]);
end

%% ===== 第9步：模糊函数分析 =====
fprintf('计算模糊函数...\n');
[ambig_func, tau_axis, nu_axis, metrics] = ambiguityFunction(s_tx, fs, c, f0);

if plot_AF
    figure('Name','模糊函数分析');
    
    % 3D模糊图
    subplot(2,2,1);
    n_tau_plot = min(500, length(tau_axis));
    n_nu_plot = min(500, length(nu_axis));
    [TAU, NU] = meshgrid(tau_axis(1:n_tau_plot)*1e6, nu_axis(1:n_nu_plot));
    surf(TAU, NU, 20*log10(ambig_func(1:n_nu_plot, 1:n_tau_plot) + eps), 'EdgeColor', 'none');
    xlabel('延迟 (μs)'); ylabel('多普勒频移 (Hz)');
    zlabel('幅度 (dB)');
    title('模糊函数 - 3D视图');
    colorbar;
    
    % 距离切片
    subplot(2,2,2);
    plot(tau_axis*1e6, 20*log10(ambig_func(round(end/2), :) + eps));
    xlabel('延迟 (μs)'); ylabel('幅度 (dB)');
    title('模糊函数 - 距离切片');
    grid on;
    
    % 多普勒切片
    subplot(2,2,3);
    plot(nu_axis, 20*log10(ambig_func(:, round(end/2)) + eps));
    xlabel('多普勒频移 (Hz)'); ylabel('幅度 (dB)');
    title('模糊函数 - 多普勒切片');
    grid on;
    
    % 性能指标
    subplot(2,2,4);
    axis off;
    text_str = sprintf(['距离分辨率: %.2f m\n', ...
                        '速度分辨率: %.2f m/s\n', ...
                        '主瓣积分宽度(τ): %.2e s\n', ...
                        '主瓣积分宽度(ν): %.2e Hz\n', ...
                        '最高旁瓣: %.2f dB'], ...
                       metrics.range_resolution, ...
                       metrics.velocity_resolution, ...
                       metrics.main_lobe_width_tau, ...
                       metrics.main_lobe_width_nu, ...
                       metrics.peak_side_lobe);
    text(0.1, 0.5, text_str, 'FontSize', 11, 'FontFamily', 'monospace');
    title('性能指标');
end

%% ===== 第10步：特征提取和目标检测 =====
fprintf('目标检测和特征提取...\n');
[detected_range, detected_velocity, detection_confidence] = featureExtract(rd_map, range_axis_plot, velocity_axis_plot, detection_threshold);

%% ===== 第11步：引信判决 =====
fprintf('引信判决...\n');
[detonation_signal, early_detonation_risk, analysis_result] = fuzeDecision(detected_range, detected_velocity, detection_confidence, R_target, R_safe, v_target);

%% ===== 第12步：输出结果 =====
fprintf('\n========== 仿真结果 ==========\n');
fprintf('真实目标参数：\n');
fprintf('  距离: %.1f m\n', R_target);
fprintf('  速度: %.1f m/s\n', v_target);
fprintf('\n检测结果：\n');
if ~isnan(detected_range)
    fprintf('  检测距离: %.1f m\n', detected_range);
    fprintf('  检测速度: %.1f m/s\n', detected_velocity);
    fprintf('  置信度: %.2f\n', detection_confidence);
    fprintf('  距离误差: %.1f m\n', detected_range - R_target);
else
    fprintf('  未检测到目标\n');
end

fprintf('\n引信状态：\n');
fprintf('  安全距离: %.1f m\n', R_safe);
fprintf('  引爆距离: %.1f m\n', R_target - R_safe);
if early_detonation_risk
    fprintf('  【警告】存在早炸风险！\n');
    fprintf('  最早检测距离: %.1f m\n', analysis_result.earliest_range);
    fprintf('  安全余度: %.1f m\n', analysis_result.safety_margin);
else
    fprintf('  引信安全\n');
end

fprintf('  起爆信号: %s\n', detonation_signal);
fprintf('==============================\n\n');

%% ===== 第13步：绘制特征提取结果 =====
if plot_RDMap
    figure('Name','目标检测结果');
    imagesc(velocity_axis_plot, range_axis_plot, 20*log10(rd_map + eps));
    colorbar;
    xlabel('速度 (m/s)'); ylabel('距离 (m)');
    title('距离-多普勒图 及 检测结果');
    axis xy;
    set(gca, 'YLim', [0 1000]);
    hold on;
    
    if ~isnan(detected_range)
        plot(detected_velocity, detected_range, 'r*', 'MarkerSize', 15, 'LineWidth', 2);
        legend('检测峰值');
    end
    
    % 标记真实目标
    plot(v_target, R_target, 'go', 'MarkerSize', 10, 'LineWidth', 2);
    plot(v_target, R_target - R_safe, 'b^', 'MarkerSize', 10, 'LineWidth', 2);
    legend('检测峰值', '真实目标', '引爆点');
end

%% ===== 保存仿真结果 =====
save('simulation_result.mat', 's_tx', 'r_echo', 'i_total', 'n_noise', 'y_rx', ...
     'y_compressed', 'rd_map', 'detected_range', 'detected_velocity', ...
     'detection_confidence', 'detonation_signal', 'early_detonation_risk');

fprintf('仿真完成！结果已保存到 simulation_result.mat\n');
