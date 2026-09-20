function [psi, dR, dz] = gs_exact(R, z, p)
%GS_EXACT  Grad-Shafranov 方程的 Solov'ev 解析解
%
%   [psi, dR, dz] = gs_exact(R, z, p)
%
% 取 Solov'ev 型剖面（p、F^2 都是 psi 的线性函数）：
%       mu0 * p(psi)  =  p0 + A * psi          =>  mu0 * dp/dpsi = A
%       F^2(psi)      =  F0^2 + F1 * psi       =>  F * dF/dpsi  = F1/2 = B
% 于是 Grad-Shafranov 方程
%       Delta* psi = -mu0*R^2*dp/dpsi - F*dF/dpsi
% 退化为右端为常系数的线性椭圆方程
%       Delta* psi = -A*R^2 - B
% 其中  Delta* = d^2/dR^2 - (1/R)*d/dR + d^2/dz^2
%
% 本函数给出的解析解：
%       psi(R,z) = A/8 * [ R0^4 - (R^2 - R0^2)^2 ] - B/2 * z^2
%
% 验证（直接代入）：
%   Delta*(R^4) = 8 R^2,  Delta*(R^2) = 0,  Delta*(1) = 0,  Delta*(z^2) = 2
%   => Delta* psi = A/8*(0 - 8 R^2) - B/2*2 = -A*R^2 - B   ✓
%
% 磁面形状：psi = const 等价于
%       (R^2 - R0^2)^2 + (4B/A) z^2 = const
% 是以 (R0, 0) 为中心的一族闭合曲线（磁轴 O 点在 (R0,0)，psi 取极大）。
%
% 输入
%   R, z : 坐标（可为数组，与本函数逐元素运算）
%   p    : 参数结构体，字段 A, B, R0, F0, F1（F0, F1 仅 q 剖面用到）
%
% 输出
%   psi  : 磁通函数
%   dR   : dpsi/dR
%   dz   : dpsi/dz

    t   = R.^2 - p.R0^2;
    psi = p.A/8 * (p.R0^4 - t.^2) - p.B/2 * z.^2;
    dR  = -p.A/2 * R .* t;
    dz  = -p.B * z;
end
