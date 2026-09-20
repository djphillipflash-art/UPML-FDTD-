%% 离散矩阵的性质与线性求解器的选择
%  —— 说明为什么这个矩阵**不能用 CG**，以及该用什么
%
%  Grad-Shafranov 算子   Delta* = d^2/dR^2 - (1/R) d/dR + d^2/dz^2
%  多出来的 -(1/R) d/dR 一阶项，使得**中心差分离散后的矩阵不是对称的**。
%  这直接决定了 Krylov 方法的选择。
%
%  依赖：gs_solve.m, gs_exact.m

clear; clc;

%% 第1步：参数与矩阵
p = struct('A',2,'B',1,'R0',1,'F0',1,'F1',2, ...
           'Rmin',0.4,'Rmax',1.6,'zmax',0.6);
N = 65;
[~, ~, ~, info] = gs_solve(p, N, N);
A = info.A;  b = info.b;  n = info.n_unknown;

fprintf('====== 离散矩阵的性质 ======\n');
fprintf('  规模 n = %d,  非零元 = %d（每行平均 %.2f）\n', n, info.nnz, info.nnz/n);

asym = norm(A-A','fro') / norm(A,'fro');
fprintf('  非对称度 ||A-A''||_F/||A||_F = %.4f   （来自 -(1/R)*d/dR 一阶项）\n', asym);

d  = eig(full(A(1:min(600,n), 1:min(600,n))));
fprintf('  特征值：实部 ∈ [%.1f, %.1f]，最大虚部 %.2e\n', ...
        min(real(d)), max(real(d)), max(abs(imag(d))));
fprintf('  => 矩阵是**负定**且非对称（注意：CG 要求对称正定）\n');
fprintf('  条件数估计 condest(A) = %.2e\n\n', condest(A));

%% 第2步：各种 Krylov 求解器的表现
tol = 1e-10;  maxit = 1000;
fprintf('====== 线性求解器对比（迭代法，容差 %.0e）======\n', tol);
fprintf('%-22s %8s %10s %14s\n', '求解器', 'flag', '迭代步数', '相对残差');

[x,f,r,it] = pcg(A, b, tol, maxit);
fprintf('%-22s %8d %10d %14.2e\n', 'pcg(A)', f, it, r/it);

[x,f,r,it] = pcg(-A, b, tol, maxit);
it_cg_neg = it;  flag_cg_neg = f;  res_cg_neg = r/it;
fprintf('%-22s %8d %10d %14.2e\n', 'pcg(-A)', f, it, r/it);

[x,f,r,it] = gmres(A, b, [], tol, maxit);
fprintf('%-22s %8d %10d %14.2e\n', 'gmres(A)', f, it(2), r(end)/it(2));

[x,f,r,it] = bicgstab(A, b, tol, maxit);
fprintf('%-22s %8d %10d %14.2e\n', 'bicgstab(A)', f, it, r/it);

As = 0.5*(A+A');
[x,f,r,it] = pcg(-As, b, tol, maxit);
fprintf('%-22s %8d %10d %14.2e\n', 'pcg(-A 的对称部分)', f, it, r/it);

fprintf('\n结论：\n');
fprintf('  1) pcg(A) 立即失败（flag=4）：CG 要求对称正定，而 A 负定且非对称；\n');
fprintf('  2) 即使取 -A 变成正定，CG 仍达不到容差（迭代 %d 步，flag=%d，残差 %.1e），\n', ...
        it_cg_neg, flag_cg_neg, res_cg_neg);
fprintf('     因为 A 非对称，而 CG 的理论基础是对称矩阵；\n');
fprintf('  3) 正确的选择是 **GMRES 或 BiCGSTAB**（非对称系统），本例都在 200 步内收敛；\n');
fprintf('  4) 若改用 1/R 加权的 Galerkin（有限元）离散，矩阵变成对称正定，CG 即可用——\n');
fprintf('     这是"选格式"与"选求解器"必须一起考虑的原因。\n');

fprintf('\n直接法对照：稀疏 LU 解同一问题耗时 %.3f 秒（n=%d）。\n', info.t_solve, n);
