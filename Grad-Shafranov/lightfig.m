function lightfig(f)
%LIGHTFIG 把图设为白底（MATLAB R2025a 起默认深色主题，导出到报告/PPT 需要白底）
%
%   lightfig(f)  f 为图窗句柄，缺省为当前图窗。

    if nargin < 1 || isempty(f), f = gcf; end

    try, theme(f, 'light'); catch, end          % R2025a+；老版本走下面的显式设置
    set(f, 'Color', 'w');

    ax = findall(f, 'Type', 'axes');
    set(ax, 'Color', 'none', 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k');
    set(findall(f, 'Type', 'text'), 'Color', 'k');   % 标题、坐标轴标签、图例文字

    lg = findall(f, 'Type', 'legend');
    if ~isempty(lg)
        set(lg, 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.4 0.4 0.4]);
    end
end
