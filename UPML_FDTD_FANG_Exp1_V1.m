%% 第1步：参数设置与网格定义（矩形散射体）
clear; clc; close all;

% ===== 物理参数（论文例1） =====
sigma = 4.0;          % PML 衰减系数
tau = 0.005;          % 时间步长
l1 = 0.1;  l2 = 0.1;  % 空间步长（正方形网格）
M = 10000;              % 总迭代步数（先设为100，用于测试，之后改为10000）

% ===== 网格索引定义 =====
I0 = 10;  J0 = 10;    % 散射体半宽（E3 网格点数）
I1 = 50;  J1 = 50;    % PML 内边界
I2 = 55;  J2 = 55;    % PML 外边界（截断边界）

% ===== 初始波参数 =====
k1 = 1;  k2 = 1;  w = 4.5;

% ===== 网格数组尺寸 =====
Nx_E = 2*I2 + 1;      % E3 的行数（对应 i 从 -I2 到 I2）
Ny_E = 2*J2 + 1;      % E3 的列数（对应 j 从 -J2 到 J2）
Nx_H1 = 2*I2 + 1;     % H1 的行数（i 从 -I2 到 I2）
Ny_H1 = 2*J2;         % H1 的列数（j 从 -J2 到 J2-1）
Nx_H2 = 2*I2;         % H2 的行数（i 从 -I2 到 I2-1）
Ny_H2 = 2*J2 + 1;     % H2 的列数（j 从 -J2 到 J2）

fprintf('网格尺寸：E3 = %d x %d, H1 = %d x %d, H2 = %d x %d\n', ...
    Nx_E, Ny_E, Nx_H1, Ny_H1, Nx_H2, Ny_H2);
%% 第2步：生成区域标志矩阵
% E3_region: 0=边界/散射体, 1=Ω1, 2=Ω2, 3=Ω3, 4=Ω4
E3_region = zeros(Nx_E, Ny_E);
for i = -I2 : I2
    for j = -J2 : J2
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        % 外边界或散射体内部（含边界）视为0（不更新）
        if (abs(i) == I2) || (abs(j) == J2) || (abs(i) <= I0 && abs(j) <= J0)
            E3_region(idx_i, idx_j) = 0;
            continue;
        end
        % 判断 PML 状态
        in_x_pml = (abs(i) > I1) && (abs(i) < I2);
        in_y_pml = (abs(j) > J1) && (abs(j) < J2);
        if ~in_x_pml && ~in_y_pml
            E3_region(idx_i, idx_j) = 1; %在Ω1中
        elseif in_x_pml && ~in_y_pml
            E3_region(idx_i, idx_j) = 2; %在Ω2中 （使用 h2, H1）
        elseif ~in_x_pml && in_y_pml
            E3_region(idx_i, idx_j) = 3; %在Ω3中 （使用 h1, H2）
        elseif in_x_pml && in_y_pml
            E3_region(idx_i, idx_j) = 4; %在Ω4中 （使用 h2, h1）
        end
    end
end

% H1_region: 0=散射体, 1=内部用E3, 12=Ω12用E3, 34=Ω34用e3
H1_region = zeros(Nx_H1, Ny_H1);
for i = -I2 : I2
    for j = -J2 : J2-1   % H1 的列索引对应物理 j+0.5
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        % 散射体内部
        if abs(i) <= I0 && abs(j) <= J0-1
            H1_region(idx_i, idx_j) = 0;
            continue;
        end
        % 只有当 x 方向在 PML 中时，H1 才可能用 e3
        if abs(i) > I1 && abs(i) < I2
            if abs(j+0.5) <= J1 % y 方向在内部
                H1_region(idx_i, idx_j) = 12; %Ω12, 用E3
            else
                H1_region(idx_i, idx_j) = 34; %Ω34, 用e3
            end
        else
            H1_region(idx_i, idx_j) = 1;      %内部，用E3
        end
    end
end

% H2_region: 0=散射体, 1=内部用E3, 13=Ω13用E3, 24=Ω24用e3
H2_region = zeros(Nx_H2, Ny_H2);
for i = -I2 : I2-1   % H2 的行索引对应物理 i+0.5
    for j = -J2 : J2
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        if abs(i) <= I0-1 && abs(j) <= J0
            H2_region(idx_i, idx_j) = 0;
            continue;
        end
        if abs(j) > J1 && abs(j) < J2       % y 方向在 PML 中
            if abs(i+0.5) <= I1             % x 方向在内部
                H2_region(idx_i, idx_j) = 13;   % Ω13
            else
                H2_region(idx_i, idx_j) = 24;   % Ω24
            end
        else
            H2_region(idx_i, idx_j) = 1;        % 内部
        end
    end
end

fprintf('区域标志矩阵生成完成。\n');
%% 第三步：初始化电磁场和辅助变量
E3 = zeros(Nx_E, Ny_E);
H1 = zeros(Nx_H1, Ny_H1);
H2 = zeros(Nx_H2, Ny_H2);

% 初始化 H1 （在 Ω1 内部，排除散射体）
for i = -I1 : I1
    for j = -J1 : J1-1    % H1 物理 j+0.5
        if abs(i) <= I0 && abs(j) <= J0-1
            continue;  % 散射体内部
        end
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        % 公式 （2.1）下方的初始条件
        H1(idx_i, idx_j) = (k1 + 0.5) * cos((k1+0.5)*pi*(i+0.5)*l1 + 0.5*pi) ...
                           * sin((k2+0.5)*pi*j*l2 + 0.5*pi) / w;
    end
end

% 初始化 H2 (在 Ω1 内部，排除散射体)
for i = -I1 : I1-1
    for j = -J1 : J1
        if abs(i) <= I0-1 && abs(j) <= J0
            continue;
        end
        idx_i = i + I2 + 1;
        idx_j = j + J2 + 1;
        H2(idx_i, idx_j) = (k1 + 0.5) * cos((k1+0.5)*pi*(i+0.5)*l1 + 0.5*pi) ...
            * sin((k2+0.5)*pi*j*l2 + 0.5*pi) / w;
    end
end

% 辅助变量初始等于原变量
h1 = H1;
h2 = H2;
e3 = E3;

% 保存上一时刻场（用于递推辅助变量）
E3_prev = E3;
H1_prev = H1;
H2_prev = H2;

% 辅助参数 α, β
alpha = 1 + sigma * tau / 2;
beta  = 1 - sigma * tau / 2;

fprintf('场初始化完成。\n');

%% 第5步：完整主循环（更新所有场和辅助变量）
fprintf('开始完整循环 (M=%d)...\n', M);
tic;

record_step = 100;
num_records = floor(M / record_step);
norms = zeros(num_records, 1);
times = zeros(num_records, 1);

for n = 1 : M
    % ---- 1. 更新 E3 (同第4步) ----
   % ---- 1. 更新 E3 ----
    for i = -I2+1 : I2-1
        for j = -J2+1 : J2-1
            idx_i = i + I2 + 1;
            idx_j = j + J2 + 1;
            region = E3_region(idx_i, idx_j);
            if region == 0
                continue;
            end
            
            % 【修正】严格对应 Yee 网格的中心差分 (物理坐标 i+0.5 和 i-0.5)
            h2_right = H2(i   + I2 + 1, idx_j);
            h2_left  = H2(i-1 + I2 + 1, idx_j);
            dH2 = (h2_right - h2_left) / l1;
            
            % 【修正】严格对应 Yee 网格的中心差分 (物理坐标 j+0.5 和 j-0.5)
            h1_up   = H1(idx_i, j   + J2 + 1);
            h1_down = H1(idx_i, j-1 + J2 + 1);
            dH1 = (h1_up - h1_down) / l2;
            
            % 根据区域选择更新变量 (同步修正辅助变量 h1, h2 的索引)
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
    % 边界强制0
    E3(1, :) = 0; E3(end, :) = 0; E3(:, 1) = 0; E3(:, end) = 0;
    for i = -I0 : I0
        for j = -J0 : J0
            if abs(i) == I0 || abs(j) == J0
                E3(i+I2+1, j+J2+1) = 0;
            end
        end
    end
    
    % ---- 2. 更新辅助变量 e3 ----
    e3 = (beta/alpha) * e3 + (1/alpha) * (E3 - E3_prev);
    
    % ---- 3. 更新 H1 ----
    for i = -I2 : I2
        for j = -J2 : J2-1
            idx_i = i + I2 + 1;
            idx_j = j + J2 + 1;
            region = H1_region(idx_i, idx_j);
            if region == 0
                continue;
            end
            % 计算 dy 差分
            if j+1 <= J2-1
                if region == 34
                    d = (e3(idx_i, j+1+J2+1) - e3(idx_i, j+J2+1)) / l2;
                else
                    d = (E3(idx_i, j+1+J2+1) - E3(idx_i, j+J2+1)) / l2;
                end
            else
                d = 0;
            end
            H1(idx_i, idx_j) = H1(idx_i, idx_j) + tau * d;
        end
    end
    
    % ---- 4. 更新 H2 ----
    for i = -I2 : I2-1
        for j = -J2 : J2
            idx_i = i + I2 + 1;
            idx_j = j + J2 + 1;
            region = H2_region(idx_i, idx_j);
            if region == 0
                continue;
            end
            if i+1 <= I2-1
                if region == 24
                    d = (e3(i+1+I2+1, idx_j) - e3(i+I2+1, idx_j)) / l1;
                else
                    d = (E3(i+1+I2+1, idx_j) - E3(i+I2+1, idx_j)) / l1;
                end
            else
                d = 0;
            end
            H2(idx_i, idx_j) = H2(idx_i, idx_j) + tau * d;
        end
    end
    
    % ---- 5. 更新辅助变量 h1, h2 ----
    h1 = (beta/alpha) * h1 + (1/alpha) * (H1 - H1_prev);
    h2 = (beta/alpha) * h2 + (1/alpha) * (H2 - H2_prev);
    
    % ---- 保存上一时刻场 ----
    E3_prev = E3;
    H1_prev = H1;
    H2_prev = H2;
    
    % ---- 记录范数 ----
    if mod(n, record_step) == 0
        sum_sq = 0;
        % H1
        for i = -I1 : I1
            for j = -J1 : J1-1
                if abs(i) <= I0 && abs(j) <= J0-1, continue; end
                idx_i = i + I2 + 1;
                idx_j = j + J2 + 1;
                sum_sq = sum_sq + H1(idx_i, idx_j)^2;
            end
        end
        % H2
        for i = -I1 : I1-1
            for j = -J1 : J1
                if abs(i) <= I0-1 && abs(j) <= J0, continue; end
                idx_i = i + I2 + 1;
                idx_j = j + J2 + 1;
                sum_sq = sum_sq + H2(idx_i, idx_j)^2;
            end
        end
        % E3
        for i = -I1 : I1
            for j = -J1 : J1
                if abs(i) <= I0 && abs(j) <= J0, continue; end
                idx_i = i + I2 + 1;
                idx_j = j + J2 + 1;
                sum_sq = sum_sq + E3(idx_i, idx_j)^2;
            end
        end
        norm_val = sqrt(l1 * l2 * sum_sq);
        norms(floor(n/record_step)) = norm_val;
        times(floor(n/record_step)) = n * tau;
    end
end

toc;
fprintf('完整循环结束。\n');

% 绘图
figure;
plot(times, norms, 'b-', 'LineWidth', 1.5);
xlabel('Time t');
ylabel('l^2(\Omega_1) norm');
title('FDTD-UPML Stability (Rectangular scatterer)');
grid on;
