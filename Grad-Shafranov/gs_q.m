function q = gs_q(rho, p, Rg, zg, psig)
%GS_Q  计算磁面 rho 上的安全因子 q
%
%   q(rho) = F/(2*pi) * ∮ dl / (R |grad psi|)
%
% 磁面的解析参数化（由 Solov'ev 解的磁面方程得到）：
%       (R^2 - R0^2)^2 + (4B/A) z^2 = rho^2
%   令  R^2 - R0^2 = rho*cos(theta),  z = (rho/2)*sqrt(A/B)*sin(theta)
%
% 输入
%   rho  : 磁面标号（0 < rho < R0^2，rho 越大磁面越大）
%   p    : 参数结构体
%   Rg,zg,psig : 数值解的网格与 psi（传空数组 [] 表示改用解析解）
%
% 输出
%   q    : 该磁面上的安全因子
%
% 说明：解析解与数值解都可以调用本函数（后者通过 interp2 插值），
%       两条 q 曲线的比较构成一层独立的验证。

    nth = 1441;
    th  = linspace(0, 2*pi, nth)';

    % 磁面参数化
    R = sqrt(p.R0^2 + rho*cos(th));
    z = (rho/2) * sqrt(p.A/p.B) * sin(th);

    if isempty(psig)
        % ---- 用解析解 ----
        [psi_f, dR_f, dz_f] = gs_exact(R, z, p);
    else
        % ---- 用数值解插值 ----
        % 数值梯度（psig 的第 1 维是 R，第 2 维是 z）
        NR = numel(Rg);  Nz = numel(zg);
        hR = Rg(2)-Rg(1); hz = zg(2)-zg(1);
        gR = zeros(NR,Nz);  gz = zeros(NR,Nz);
        gR(2:end-1,:) = (psig(3:end,:) - psig(1:end-2,:)) / (2*hR);
        gz(:,2:end-1) = (psig(:,3:end) - psig(:,1:end-2)) / (2*hz);

        % interp2(X,Y,V,...) 要求 V 的尺寸为 (numel(Y), numel(X))
        % 这里声明 X=zg、Y=Rg，故 V 应为 numel(Rg) x numel(zg) = psig 本身
        F = @(M) interp2(zg, Rg, M, z, R, 'spline');
        psi_f = F(psig);
        dR_f  = F(gR);
        dz_f  = F(gz);
    end

    gradpsi = sqrt(dR_f.^2 + dz_f.^2);

    % 弧长元
    dR = gradient(R, th);
    dz = gradient(z, th);
    dl = sqrt(dR.^2 + dz.^2);

    % F(psi) = sqrt(F0^2 + F1*psi)
    Fpsi = sqrt(p.F0^2 + p.F1 * psi_f);

    integrand = dl ./ (R .* gradpsi);
    q = mean(Fpsi) / (2*pi) * trapz(th, integrand);
end
