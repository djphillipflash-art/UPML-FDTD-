%% 二维 FDTD-UPML 数值稳定性：能量范数随时间的变化
%  —— 复现论文 Fang & Ying (2009) §6 例 1（矩形散射体）
%
%  计算域：[-a-d, a+d] × [-a-d, a+d]
%      Ω1 = |x1| <= a 且 |x2| <= a          物理区
%      PML = a < |x1| < a+d 或 a < |x2| < a+d
%      外边界 |x1| = a+d 或 |x2| = a+d 上 E3 = 0
%      散射体 |x1| <= I0*l 且 |x2| <= J0*l 上 E3 = 0（矩形理想导体）
%
%  场量布局（Yee 交错网格）：
%      E3 在整点  (i, j)      、整时刻  n
%      H1 在半点  (i, j+1/2)  、半时刻 n+1/2
%      H2 在半点  (i+1/2, j)  、半时刻 n+1/2
%
%  UPML 辅助变量（论文 (3.3) 的 b/g 变量）：
%      alpha = 1 + sigma*tau/2,  beta = 1 - sigma*tau/2
%      递推：e3 = (beta/alpha)*e3 + (1/alpha)*(E3 - E3_prev)，h1、h2 同理
%
%  区域切换（论文 (3.4)—(3.11)）：
%      H1 在 x2 位于 PML 时用 e3，否则用 E3        （Ω34 / Ω12）
%      H2 在 x1 位于 PML 时用 e3，否则用 E3        （Ω24 / Ω13）
%      E3 在 x1 位于 PML 时用 h2，在 x2 位于 PML 时用 h1
%
%  初始条件（论文 §6 例 1）：
%      H1^{1/2}(i,j+1/2)
%          = (k2+0.5)*sin((k1+0.5)*pi*i*l1 + 0.5*pi)*cos((k2+0.5)*pi*j*l2 + 0.5*pi)/w
%      H2^{1/2}(i+1/2,j)
%          = (k1+0.5)*cos((k1+0.5)*pi*i*l1 + 0.5*pi)*sin((k2+0.5)*pi*j*l2 + 0.5*pi)/w
%      仅在 Ω1 内非零；E3^0 = 0，E3^1 由差分方程给出（即第一步迭代）

clear; clc; close all;

%% 第1步：参数设置与网格定义（矩形散射体）
sigma = 4.0;         % PML 衰减系数
tau   = 0.005;       % 时间步长
l1 = 0.1;  l2 = 0.1; % 空间步长（正方形网格）
M = 10000;           % 总迭代步数（T = M*tau = 50）

I0 = 10;  J0 = 10;   % 散射体半宽（E3 网格点数）
I1 = 50;  J1 = 50;   % PML 内边界
I2 = 55;  J2 = 55;   % PML 外边界（截断边界）

k1 = 1;  k2 = 1;  w = 4.5;   % 初始波参数

record_step = 10;    % 每隔多少步记录一次范数

% 网格数组尺寸
Nx_E = 2*I2 + 1;     % E3 的行数（i 从 -I2 到 I2）
Ny_E = 2*J2 + 1;     % E3 的列数（j 从 -J2 到 J2）
Nx_H1 = 2*I2 + 1;    % H1 的行数
Ny_H1 = 2*J2;        % H1 的列数（j 从 -J2 到 J2-1）
Nx_H2 = 2*I2;        % H2 的行数（i 从 -I2 到 I2-1）
Ny_H2 = 2*J2 + 1;    % H2 的列数

fprintf('二维 FDTD-UPML：sigma = %g, tau = %g, l1 = l2 = %g, M = %d\n', sigma, tau, l1, M);
fprintf('网格尺寸：E3 = %d x %d, H1 = %d x %d, H2 = %d x %d\n', ...
    Nx_E, Ny_E, Nx_H1, Ny_H1, Nx_H2, Ny_H2);
fprintf('物理区 |x| <= %.1f，PML 宽 %.1f（%d 层）\n', I1*l1, (I2-I1)*l1, I2-I1);
fprintf('库朗数 tau/l = %.4f（二维 CFL 上限 %.4f）\n', tau/l1, 1/sqrt(2));

%% 第2步：生成区域标志矩阵
% E3_region：0=边界/散射体, 1=Ω1, 2=Ω2, 3=Ω3, 4=Ω4
E3_region = zeros(Nx_E, Ny_E);
for i = -I2 : I2
    for j = -J2 : J2
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        % 外边界或散射体内部（含边界）视为 0（不更新）
        if (abs(i) == I2) || (abs(j) == J2) || (abs(i) <= I0 && abs(j) <= J0)
            E3_region(idx_i, idx_j) = 0;
            continue;
        end
        in_x_pml = (abs(i) > I1) && (abs(i) < I2);
        in_y_pml = (abs(j) > J1) && (abs(j) < J2);
        if ~in_x_pml && ~in_y_pml
            E3_region(idx_i, idx_j) = 1;   % Ω1
        elseif in_x_pml && ~in_y_pml
            E3_region(idx_i, idx_j) = 2;   % Ω2：x1 在 PML，用 h2
        elseif ~in_x_pml && in_y_pml
            E3_region(idx_i, idx_j) = 3;   % Ω3：x2 在 PML，用 h1
        else
            E3_region(idx_i, idx_j) = 4;   % Ω4：两个方向都在 PML
        end
    end
end

% H1_region：0=散射体, 12=Ω12（用 E3）, 34=Ω34（用 e3）
%   注意：只由 x2 是否在 PML 决定，与 x1 无关（论文 Ω34 = {x2 在 PML}）
H1_region = zeros(Nx_H1, Ny_H1);
for i = -I2 : I2
    for j = -J2 : J2-1
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        if abs(i) <= I0 && abs(j) <= J0-1
            H1_region(idx_i, idx_j) = 0;   % 散射体内部
            continue;
        end
        if abs(j+0.5) > J1 && abs(j+0.5) < J2
            H1_region(idx_i, idx_j) = 34;  % x2 在 PML，用 e3
        else
            H1_region(idx_i, idx_j) = 12;  % x2 在 Ω1，用 E3
        end
    end
end

% H2_region：0=散射体, 13=Ω13（用 E3）, 24=Ω24（用 e3）
%   注意：只由 x1 是否在 PML 决定，与 x2 无关（论文 Ω24 = {x1 在 PML}）
H2_region = zeros(Nx_H2, Ny_H2);
for i = -I2 : I2-1
    for j = -J2 : J2
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        if abs(i) <= I0-1 && abs(j) <= J0
            H2_region(idx_i, idx_j) = 0;   % 散射体内部
            continue;
        end
        if abs(i+0.5) > I1 && abs(i+0.5) < I2
            H2_region(idx_i, idx_j) = 24;  % x1 在 PML，用 e3
        else
            H2_region(idx_i, idx_j) = 13;  % x1 在 Ω1，用 E3
        end
    end
end

fprintf('区域标志矩阵生成完成。\n');

%% 第3步：初始化电磁场与辅助变量
E3 = zeros(Nx_E, Ny_E);
H1 = zeros(Nx_H1, Ny_H1);
H2 = zeros(Nx_H2, Ny_H2);

% 初始化 H1（仅在 Ω1 内，排除散射体）
for i = -I1 : I1
    for j = -J1 : J1-1
        if abs(i) <= I0 && abs(j) <= J0-1
            continue;
        end
        H1(i+I2+1, j+J2+1) = (k2 + 0.5) * sin((k1+0.5)*pi*i*l1 + 0.5*pi) ...
                           * cos((k2+0.5)*pi*j*l2 + 0.5*pi) / w;
    end
end

% 初始化 H2（仅在 Ω1 内，排除散射体）
for i = -I1 : I1-1
    for j = -J1 : J1
        if abs(i) <= I0-1 && abs(j) <= J0
            continue;
        end
        H2(i+I2+1, j+J2+1) = (k1 + 0.5) * cos((k1+0.5)*pi*i*l1 + 0.5*pi) ...
                           * sin((k2+0.5)*pi*j*l2 + 0.5*pi) / w;
    end
end

% 辅助变量初始等于原变量
e3 = E3;
h1 = H1;
h2 = H2;

% 保存上一时刻场（用于递推辅助变量）
E3_prev = E3;
H1_prev = H1;
H2_prev = H2;

% 辅助参数
alpha = 1 + sigma * tau / 2;
beta  = 1 - sigma * tau / 2;

fprintf('场初始化完成。\n');

%% 第4步：主循环
num_records = floor(M / record_step);
norms = zeros(num_records, 1);
times = zeros(num_records, 1);

fprintf('开始主循环 (M=%d)...\n', M);
tic;

for n = 1 : M

    % ---- 4.1 更新 E3 ----
    for i = -I2+1 : I2-1
        for j = -J2+1 : J2-1
            idx_i = i + I2 + 1;
            idx_j = j + J2 + 1;
            region = E3_region(idx_i, idx_j);
            if region == 0
                continue;
            end

            % x1 方向的差分：H2 在 (i±1/2, j)
            dH2 = (H2(i+I2+1, idx_j) - H2(i-1+I2+1, idx_j)) / l1;
            % x2 方向的差分：H1 在 (i, j±1/2)
            dH1 = (H1(idx_i, j+J2+1) - H1(idx_i, j-1+J2+1)) / l2;

            % 根据区域选择用原变量还是辅助变量（论文 (3.8)—(3.11)）
            switch region
                case 1
                    d2 = dH2;  d1 = dH1;
                case 2
                    d2 = (h2(i+I2+1, idx_j) - h2(i-1+I2+1, idx_j)) / l1;
                    d1 = dH1;
                case 3
                    d2 = dH2;
                    d1 = (h1(idx_i, j+J2+1) - h1(idx_i, j-1+J2+1)) / l2;
                case 4
                    d2 = (h2(i+I2+1, idx_j) - h2(i-1+I2+1, idx_j)) / l1;
                    d1 = (h1(idx_i, j+J2+1) - h1(idx_i, j-1+J2+1)) / l2;
            end
            E3(idx_i, idx_j) = E3(idx_i, idx_j) + tau * (d2 + d1);
        end
    end

    % 外边界强制为 0
    E3(1, :) = 0;  E3(end, :) = 0;  E3(:, 1) = 0;  E3(:, end) = 0;

    % 散射体边界强制为 0
    for i = -I0 : I0
        for j = -J0 : J0
            if abs(i) == I0 || abs(j) == J0
                E3(i+I2+1, j+J2+1) = 0;
            end
        end
    end

    % ---- 4.2 更新辅助变量 e3 ----
    e3 = (beta/alpha) * e3 + (1/alpha) * (E3 - E3_prev);

    % ---- 4.3 更新 H1 ----
    for i = -I2 : I2
        for j = -J2 : J2-1
            idx_i = i + I2 + 1;
            idx_j = j + J2 + 1;
            region = H1_region(idx_i, idx_j);
            if region == 0
                continue;
            end
            if region == 34
                d = (e3(idx_i, j+1+J2+1) - e3(idx_i, j+J2+1)) / l2;
            else
                d = (E3(idx_i, j+1+J2+1) - E3(idx_i, j+J2+1)) / l2;
            end
            H1(idx_i, idx_j) = H1(idx_i, idx_j) + tau * d;
        end
    end

    % ---- 4.4 更新 H2 ----
    for i = -I2 : I2-1
        for j = -J2 : J2
            idx_i = i + I2 + 1;
            idx_j = j + J2 + 1;
            region = H2_region(idx_i, idx_j);
            if region == 0
                continue;
            end
            if region == 24
                d = (e3(i+1+I2+1, idx_j) - e3(i+I2+1, idx_j)) / l1;
            else
                d = (E3(i+1+I2+1, idx_j) - E3(i+I2+1, idx_j)) / l1;
            end
            H2(idx_i, idx_j) = H2(idx_i, idx_j) + tau * d;
        end
    end

    % ---- 4.5 更新辅助变量 h1, h2 ----
    h1 = (beta/alpha) * h1 + (1/alpha) * (H1 - H1_prev);
    h2 = (beta/alpha) * h2 + (1/alpha) * (H2 - H2_prev);

    % ---- 4.6 保存上一时刻场 ----
    E3_prev = E3;
    H1_prev = H1;
    H2_prev = H2;

    % ---- 4.7 记录 Ω1 内的能量范数 ----
    if mod(n, record_step) == 0
        sum_sq = 0;
        for i = -I1 : I1
            for j = -J1 : J1
                if abs(i) <= I0 && abs(j) <= J0
                    continue;
                end
                sum_sq = sum_sq + E3(i+I2+1, j+J2+1)^2;
            end
        end
        for i = -I1 : I1
            for j = -J1 : J1-1
                if abs(i) <= I0 && abs(j) <= J0-1
                    continue;
                end
                sum_sq = sum_sq + H1(i+I2+1, j+J2+1)^2;
            end
        end
        for i = -I1 : I1-1
            for j = -J1 : J1
                if abs(i) <= I0-1 && abs(j) <= J0
                    continue;
                end
                sum_sq = sum_sq + H2(i+I2+1, j+J2+1)^2;
            end
        end
        norms(floor(n/record_step)) = sqrt(l1 * l2 * sum_sq);
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
title(sprintf('2D FDTD-UPML: energy norm vs time  (\\sigma = %g, S = %.3f)', sigma, tau/l1));
% 设为白底（MATLAB R2025a 起默认深色主题，导出到报告/PPT 需要白底）
set(gcf, 'Color', 'w');
set(gca, 'Color', 'none', 'XColor', 'k', 'YColor', 'k', ...
    'GridColor', [0.15 0.15 0.15], 'GridAlpha', 0.15);

exportgraphics(gcf, 'norm_2d.png', 'Resolution', 150);
fprintf('\n图已保存为 norm_2d.png\n');

% 说明：范数在 t ≈ 14 后停在初值的 18% 左右，这不是数值不稳定，而是论文给的初值
%       不满足 div(B) = 0（实测 max|div B| = 2.88，而 max|H| = 0.33），
%       其中非辐射的那部分既不传播也不被 PML 吸收，会一直留在 Ω1 内。
%       把 H2 的初值反号（k1 = k2 时即满足 div(B) = 0）即可看到范数继续衰减。
