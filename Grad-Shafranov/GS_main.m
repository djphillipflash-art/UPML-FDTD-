%% Grad-Shafranov 平衡求解器：二阶收敛性验证与物理量提取
%  —— 用 Solov'ev 解析解做代码验证（code verification）
%
%  方程： Delta* psi = -A*R^2 - B
%          Delta* = d^2/dR^2 - (1/R)*d/dR + d^2/dz^2
%  解析解： psi = A/8*[R0^4 - (R^2-R0^2)^2] - B/2*z^2
%          （对应 Solov'ev 型剖面 mu0*dp/dpsi = A，F*dF/dpsi = B）
%
%  流程： 第1步 参数
%         第2步 基准网格求解
%         第3步 网格收敛性（验证二阶）
%         第4步 Richardson 外推与 GCI
%         第5步 物理量：磁轴位置、q 剖面
%         第6步 绘图与输出
%
%  依赖：gs_exact.m, gs_solve.m, gs_q.m

clear; clc; close all;

% 把本文件所在目录加入路径，保证能调用同目录的函数
thisdir = fileparts(mfilename('fullpath'));
if ~isempty(thisdir), addpath(thisdir); end

%% 第1步：参数设置
% Solov'ev 剖面参数（归一化单位）
A  = 2.0;      % mu0 * dp/dpsi
B  = 1.0;      % F * dF/dpsi
R0 = 1.0;      % 磁轴大半径
F0 = 1.0;      % F^2 = F0^2 + F1*psi 中的常数项
F1 = 2*B;      % 与 B 一致

% 计算域（关于 R0 与 z=0 对称，便于把磁轴放在网格节点上）
Rmin = R0 - 0.6;   Rmax = R0 + 0.6;
zmax = 0.6;

p = struct('A',A,'B',B,'R0',R0,'F0',F0,'F1',F1, ...
           'Rmin',Rmin,'Rmax',Rmax,'zmax',zmax);

fprintf('======== Grad-Shafranov 平衡求解器 ========\n');
fprintf('  A = mu0*dp/dpsi = %.2f,  B = F*dF/dpsi = %.2f,  R0 = %.2f\n', A, B, R0);
fprintf('  计算域： R ∈ [%.2f, %.2f],  z ∈ [%.2f, %.2f]\n', Rmin, Rmax, -zmax, zmax);
fprintf('  磁面方程： (R^2-R0^2)^2 + (4B/A) z^2 = rho^2，  rho 越大磁面越大\n');
fprintf('===========================================\n\n');

%% 第2步：基准网格求解
N0 = 129;
[psi0, R0g, z0g, info0] = gs_solve(p, N0, N0);
fprintf('基准网格 %dx%d： 未知量 %d 个， 非零元 %d 个， 求解耗时 %.3f 秒\n', ...
        N0, N0, info0.n_unknown, info0.nnz, info0.t_solve);

% 数值解与解析解的场误差
[PSIex0, ~, ~] = gs_exact(R0g*ones(1,N0), ones(N0,1)*z0g', p);
err0 = psi0 - PSIex0;
fprintf('  基准网格上： max|误差| = %.3e,  RMS 误差 = %.3e\n', ...
        max(abs(err0(:))), sqrt(mean(err0(:).^2)));
fprintf('  psi 量级： max|psi| = %.4f\n\n', max(abs(PSIex0(:))));

%% 第3步：网格收敛性（验证二阶）
Nlist = [17 33 65 129 257];
ncase = numel(Nlist);
h  = zeros(ncase,1);
e2 = zeros(ncase,1);    % RMS 误差
einf = zeros(ncase,1);  % 最大模误差
psi_axis = zeros(ncase,1);
R_axis   = zeros(ncase,1);
tsolve   = zeros(ncase,1);

fprintf('---- 网格收敛性 ----\n');
fprintf('%6s %10s %14s %14s %14s %10s\n', 'N', 'h', 'RMS 误差', 'max 误差', 'psi(轴)', '求解(s)');
for k = 1:ncase
    N = Nlist(k);
    [psik, Rk, zk, infok] = gs_solve(p, N, N);
    [PSIexk, ~, ~] = gs_exact(Rk*ones(1,N), ones(N,1)*zk', p);
    ek = psik - PSIexk;

    h(k)    = infok.hR;
    e2(k)   = sqrt(mean(ek(:).^2));
    einf(k) = max(abs(ek(:)));
    tsolve(k) = infok.t_solve;

    % 磁轴：数值解的最大值位置（网格对称，R0 与 z=0 恰为节点）
    [mx, idx] = max(psik(:));
    [ia, ja] = ind2sub(size(psik), idx);
    psi_axis(k) = mx;
    R_axis(k)   = Rk(ia);

    fprintf('%6d %10.5f %14.4e %14.4e %14.8f %10.3f\n', ...
            N, h(k), e2(k), einf(k), psi_axis(k), tsolve(k));
end

% 收敛阶 p = log(e(h1)/e(h2)) / log(h1/h2)
fprintf('\n---- 收敛阶（相邻两次加密）----\n');
fprintf('%14s %14s %10s\n', '区间', 'order(RMS)', 'order(max)');
order2 = zeros(ncase-1,1);  orderinf = zeros(ncase-1,1);
for k = 1:ncase-1
    order2(k)   = log(e2(k)/e2(k+1))   / log(h(k)/h(k+1));
    orderinf(k) = log(einf(k)/einf(k+1)) / log(h(k)/h(k+1));
    fprintf('%6d ->%6d %14.4f %10.4f\n', Nlist(k), Nlist(k+1), order2(k), orderinf(k));
end
fprintf('理论值：二阶格式，收敛阶 = 2\n');

fprintf('\n---- 磁轴位置与轴上 psi ----\n');
[~, ~, ~] = gs_exact(p.R0, 0, p);
psi_axis_exact = gs_exact(p.R0, 0, p);
fprintf('  解析解： R_axis = %.6f,  psi_axis = %.8f\n', p.R0, psi_axis_exact);
for k = 1:ncase
    fprintf('  N=%4d： R_axis = %.6f (误差 %+.2e),  psi_axis = %.8f (误差 %+.2e)\n', ...
            Nlist(k), R_axis(k), R_axis(k)-p.R0, psi_axis(k), psi_axis(k)-psi_axis_exact);
end

%% 第4步：Richardson 外推与 GCI
% 目标量：轴上 psi（加密比 r = 2，理论阶 pth = 2）
r = 2;  Fs = 1.25;  pth = 2;
f1 = psi_axis(end-1);   f2 = psi_axis(end);     % 细网格 f2、次细网格 f1
f_extrap = f2 + (f2 - f1)/(r^pth - 1);          % Richardson 外推
e_rel    = abs((f1 - f2)/f2);
GCI      = Fs * e_rel / (r^pth - 1);

fprintf('\n---- Richardson 外推与 GCI（目标量：轴上 psi）----\n');
fprintf('  次细网格 N=%d： %.10f\n', Nlist(end-1), f1);
fprintf('  细网格   N=%d： %.10f\n', Nlist(end),   f2);
fprintf('  外推值        ： %.10f\n', f_extrap);
fprintf('  解析值        ： %.10f\n', psi_axis_exact);
fprintf('  外推值相对误差： %.3e\n', abs(f_extrap-psi_axis_exact)/abs(psi_axis_exact));
fprintf('  GCI(细网格)   ： %.4f %%\n', 100*GCI);

%% 第5步：物理量——q 剖面
% 磁面标号 rho 的取值范围受计算域限制（磁面必须落在矩形域内）
Rmin_ = Rmin;  Rmax_ = Rmax;  zmax_ = zmax;
rho_lim = min([ R0^2 - Rmin_^2, Rmax_^2 - R0^2, (2*zmax_/sqrt(A/B))^2 ]);
rho_list = linspace(0.15, 0.85, 12) * sqrt(rho_lim);   % 避开磁轴附近（q→∞）

% 插值朝向自检：把解析解放到网格上再插值回来，应与解析解一致
[PSIchk, ~, ~] = gs_exact(R0g*ones(1,N0), ones(N0,1)*z0g', p);
Rq = 0.85;  zq = 0.25;
v_int = interp2(z0g, R0g, PSIchk, zq, Rq, 'spline');
v_ana = gs_exact(Rq, zq, p);
fprintf('\n插值自检： interp2 结果 = %.8f,  解析解 = %.8f,  差 = %.2e\n', ...
        v_int, v_ana, abs(v_int-v_ana));
if abs(v_int - v_ana) > 1e-4
    warning('插值朝向可能不对，请检查 interp2 的参数顺序');
end

fprintf('\n---- q 剖面（数值解 vs 解析解）----\n');
fprintf('%10s %14s %14s %12s\n', 'rho', 'q (数值)', 'q (解析)', '相对差');
q_num = zeros(size(rho_list));  q_exa = zeros(size(rho_list));
for k = 1:numel(rho_list)
    rr = rho_list(k);
    q_num(k) = gs_q(rr, p, R0g, z0g, psi0);   % 数值解（基准网格）
    q_exa(k) = gs_q(rr, p, [], [], []);       % 解析解
    fprintf('%10.4f %14.6f %14.6f %11.3f%%\n', ...
            rr, q_num(k), q_exa(k), 100*abs(q_num(k)-q_exa(k))/abs(q_exa(k)));
end

%% 第6步：绘图与输出
% ---- 图 1：磁面 ----
f1fig = figure('Visible','off','Position',[100 100 820 620]);
psip = psi0;  Rp = R0g;  zp = z0g;      % 复用第 2 步的基准网格解
[PSIexp, ~, ~] = gs_exact(Rp*ones(1,N0), ones(N0,1)*zp', p);
lev = linspace(0.05, 0.24, 14);
contour(Rp, zp, psip.',   lev, 'b-',  'LineWidth', 1.8); hold on;   % 数值解（粗蓝实线）
contour(Rp, zp, PSIexp.', lev, 'k--', 'LineWidth', 1.0);            % 解析解（细黑虚线，画在上层）
plot(p.R0, 0, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
grid on; box on;
xlabel('R'); ylabel('z');
title('Grad-Shafranov 磁面：数值解（蓝实线）vs Solov''ev 解析解（黑虚线）');
legend({'数值解', '解析解', '磁轴'}, 'Location', 'eastoutside');
lightfig(gcf);
exportgraphics(f1fig, 'fig_gs_flux.png', 'Resolution', 150);

% ---- 图 2：收敛阶 ----
f2fig = figure('Visible','off','Position',[100 100 820 560]);
loglog(h, e2,   'o-', 'LineWidth', 1.8, 'MarkerSize', 6); hold on;
loglog(h, einf, 's-', 'LineWidth', 1.5, 'MarkerSize', 5);
href = h(end)*[1 8];
loglog(href, e2(end)*(href/h(end)).^2, 'k--', 'LineWidth', 1.2);
grid on; box on;
set(gca,'XDir','reverse');
xlabel('网格步长 h'); ylabel('误差');
title(sprintf('网格收敛性：拟合阶 %.2f（RMS）/ %.2f（max），理论二阶', ...
      order2(end), orderinf(end)));
legend({'RMS 误差','max 误差','斜率 2 参考线'}, 'Location', 'northwest');
lightfig(gcf);
exportgraphics(f2fig, 'fig_gs_order.png', 'Resolution', 150);

% ---- 图 3：q 剖面 ----
f3fig = figure('Visible','off','Position',[100 100 820 560]);
plot(rho_list, q_exa, 'k--', 'LineWidth', 1.8); hold on;
plot(rho_list, q_num, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 5);
grid on; box on;
xlabel('磁面标号 \rho'); ylabel('安全因子 q');
title('安全因子剖面：数值解 vs 解析解');
legend({'解析解','数值解'}, 'Location', 'northwest');
lightfig(gcf);
exportgraphics(f3fig, 'fig_gs_q.png', 'Resolution', 150);

% ---- 数据输出 ----
T1 = table(Nlist(:), h, e2, einf, [order2; NaN], psi_axis, R_axis, tsolve, ...
     'VariableNames', {'N','h','err_RMS','err_max','order','psi_axis','R_axis','t_solve'});
writetable(T1, 'out_gs_convergence.csv', 'Encoding', 'UTF-8');

T2 = table(rho_list(:), q_num(:), q_exa(:), ...
     'VariableNames', {'rho','q_numerical','q_exact'});
writetable(T2, 'out_gs_q.csv', 'Encoding', 'UTF-8');

fprintf('\n已输出：fig_gs_flux.png, fig_gs_order.png, fig_gs_q.png, ');
fprintf('out_gs_convergence.csv, out_gs_q.csv\n');
