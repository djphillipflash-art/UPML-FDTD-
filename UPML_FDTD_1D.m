%% 一维 FDTD-UPML 数值稳定性：能量范数随时间的变化
%  —— 对应论文 Fang & Ying (2009) §6 例 1 的一维版本
%
%  计算域：x ∈ [-(a+d), a+d]
%      |x| <= a        物理区 Ω1（a = I1*l）
%      a < |x| < a+d   PML 吸收层（d = (I2-I1)*l）
%      |x| = a+d       截断边界，E = 0
%
%  场量布局（Yee 交错网格）：
%      E  在整点   i、  整时刻  n
%      H  在半点 i+1/2、半时刻 n+1/2
%
%  UPML 辅助变量（论文里的 b/g 变量）：
%      alpha = 1 + sigma*tau/2,  beta = 1 - sigma*tau/2
%      e = (beta/alpha)*e + (1/alpha)*(E - E_prev)
%      h = (beta/alpha)*h + (1/alpha)*(H - H_prev)
%      区域切换：PML 内 E 用 h 差分、H 用 e 差分；Ω1 内直接用对方原变量
%
%  初始条件：H^{1/2}(i+1/2) = (k1+0.5)*sin((k1+0.5)*pi*i*l + 0.5*pi)/w  （仅在 Ω1 内）
%            E^0 = 0，E^1 由差分方程给出（即第一步迭代）

clear; clc; close all;

%% 第1步：参数设置
sigma = 4.0;         % PML 衰减系数
tau   = 0.005;       % 时间步长
l     = 0.1;         % 空间步长
M     = 10000;       % 总迭代步数（T = M*tau = 50）

I1 = 50;             % PML 内边界（物理区半宽 a = I1*l = 5.0）
I2 = 55;             % PML 外边界（a+d = I2*l = 5.5）

k1 = 1;              % 初始波参数
w  = 4.5;

record_step = 10;    % 每隔多少步记录一次范数

fprintf('一维 FDTD-UPML：sigma = %g, tau = %g, l = %g, M = %d\n', sigma, tau, l, M);
fprintf('物理区 |x| <= %.1f，PML 宽 %.1f（%d 层），外边界 |x| = %.1f\n', ...
        I1*l, (I2-I1)*l, I2-I1, I2*l);
fprintf('库朗数 S = tau/l = %.4f（一维 CFL 上限为 1）\n', tau/l);

%% 第2步：网格与区域标志
NE = 2*I2 + 1;       % E 的点数（i = -I2 : I2）
NH = 2*I2;           % H 的点数（i = -I2 : I2-1，物理坐标 i+1/2）

% E_region：0 = 外边界（不更新），1 = Ω1（用 H 差分），2 = PML（用 h 差分）
E_region = zeros(NE, 1);
for i = -I2 : I2
    if abs(i) == I2
        continue;                        % 外边界
    end
    if abs(i) > I1
        E_region(i+I2+1) = 2;            % PML
    else
        E_region(i+I2+1) = 1;            % Ω1
    end
end

% H_region：1 = Ω1（用 E 差分），2 = PML（用 e 差分）
H_region = zeros(NH, 1);
for i = -I2 : I2-1
    if abs(i+0.5) > I1
        H_region(i+I2+1) = 2;            % PML
    else
        H_region(i+I2+1) = 1;            % Ω1
    end
end

fprintf('区域标志完成：E 中 %d 点在 Ω1、%d 点在 PML\n', ...
        sum(E_region==1), sum(E_region==2));

%% 第3步：初始化电磁场与辅助变量
E = zeros(NE, 1);
H = zeros(NH, 1);

for i = -I2 : I2-1
    if abs(i+0.5) <= I1
        H(i+I2+1) = (k1 + 0.5) * sin((k1+0.5)*pi*i*l + 0.5*pi) / w;
    end
end

e = E;               % 辅助变量初值等于原变量
h = H;
E_prev = E;          % 保存上一时刻场，用于递推辅助变量
H_prev = H;

alpha = 1 + sigma*tau/2;
beta  = 1 - sigma*tau/2;

%% 第4步：主循环
num_records = floor(M / record_step);
norms = zeros(num_records, 1);
times = zeros(num_records, 1);

fprintf('开始主循环 (M=%d)...\n', M);
tic;

for n = 1 : M

    % ---- 4.1 更新 E ----
    for i = -I2+1 : I2-1
        idx = i + I2 + 1;
        region = E_region(idx);
        if region == 0
            continue;
        end
        if region == 2
            dH = (h(idx) - h(idx-1)) / l;        % PML：用辅助变量 h
        else
            dH = (H(idx) - H(idx-1)) / l;        % Ω1：用 H
        end
        E(idx) = E(idx) + tau * dH;
    end
    E(1) = 0;  E(end) = 0;                       % 外边界强制为 0

    % ---- 4.2 更新辅助变量 e ----
    e = (beta/alpha)*e + (1/alpha)*(E - E_prev);

    % ---- 4.3 更新 H ----
    for i = -I2 : I2-1
        idx = i + I2 + 1;
        region = H_region(idx);
        if region == 2
            dE = (e(idx+1) - e(idx)) / l;        % PML：用辅助变量 e
        else
            dE = (E(idx+1) - E(idx)) / l;        % Ω1：用 E
        end
        H(idx) = H(idx) + tau * dE;
    end

    % ---- 4.4 更新辅助变量 h ----
    h = (beta/alpha)*h + (1/alpha)*(H - H_prev);

    % ---- 4.5 保存上一时刻场 ----
    E_prev = E;
    H_prev = H;

    % ---- 4.6 记录 Ω1 内的能量范数 ----
    if mod(n, record_step) == 0
        sum_sq = 0;
        for i = -I1 : I1
            sum_sq = sum_sq + E(i+I2+1)^2;
        end
        for i = -I1 : I1-1
            sum_sq = sum_sq + H(i+I2+1)^2;
        end
        norms(floor(n/record_step)) = sqrt(l * sum_sq);
        times(floor(n/record_step)) = n * tau;
    end
end

toc;
fprintf('主循环结束。\n');

%% 第5步：输出结果与绘图
fprintf('\n能量范数  l^2(Ω1)：\n');
fprintf('  初值 = %.4f\n', norms(1));
fprintf('  末值 = %.4e\n', norms(end));
fprintf('  末值/初值 = %.3e\n', norms(end)/norms(1));

figure;
semilogy(times, norms/norms(1), 'o-', 'LineWidth', 1.5, 'MarkerSize', 3);
grid on; box on;
ylim([1e-3 1.5]);          % 固定纵轴范围
xlabel('Time  t');
ylabel('Normalized energy norm  ||V||_{l^2(\Omega_1)}');
title(sprintf('1D FDTD-UPML: energy norm vs time  (\\sigma = %g, S = %.3f)', sigma, tau/l));
% 设为白底
set(gcf, 'Color', 'w');
set(gca, 'Color', 'none', 'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.15);

exportgraphics(gcf, 'norm_1d.png', 'Resolution', 150);
fprintf('\n图已保存为 norm_1d.png\n');
