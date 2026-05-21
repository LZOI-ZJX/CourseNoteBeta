% 拓扑正弦曲线绘制脚本
% 定义函数: y = sin(1/x)

% 1. 设置非均匀采样点
% 由于函数在 x 趋近于 0 时震荡无限快，使用 logspace 确保原点附近有足够的点
% 采样 100,000 个点以获得平滑的视觉效果
x = logspace(-4, 0, 100000); 
y = sin(1./x);

% 2. 绘制曲线
figure('Color', 'w');
plot(x, y, 'LineWidth', 0.5, 'Color', [0 0.4470 0.7410]);
hold on;

% 3. 绘制极限线段 {0} x [-1, 1]
% 在拓扑学中，这条线段是曲线闭包的一部分，展示了连通但不道路连通的特性
line([0 0], [-1 1], 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'Limit Segment');

% 4. 图形修饰
title('Topologist''s Sine Curve: $y = \sin(1/x)$', 'Interpreter', 'latex', 'FontSize', 14);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
grid on;
set(gca, 'XScale', 'linear'); % 如果想看极细微的震荡，可改为 'log'
xlim([-0.05 1]);
ylim([-1.2 1.2]);
legend('y = sin(1/x)', 'Limit Segment ({0} x [-1,1])', 'Location', 'northeast');

hold off;