%% 积分计算：Euler 方法与改进 Euler 方法
clear; clc;

% 1. 参数设置
h = 0.1;                % 步长
x_end = 1.0;            % 目标点
x = 0:h:x_end;          % 离散点
n = length(x);

% 初始化结果向量
y_euler = zeros(1, n);      % 标准 Euler
y_imp_euler = zeros(1, n);  % 改进 Euler

% 初值条件 I(0) = 0
y_euler(1) = 0;
y_imp_euler(1) = 0;

% 定义导数函数 f(x, y) = e^(x^2)
% 注意：此函数的导数仅与 x 有关
f = @(x) exp(x^2);

% 2. 迭代计算
for i = 1:(n-1)
    xi = x(i);
    xi_next = x(i+1);
    
    % --- 标准 Euler 方法 ---
    y_euler(i+1) = y_euler(i) + h * f(xi);
    
    % --- 改进 Euler 方法 ---
    % 预估步 (Predictor)
    yp = y_imp_euler(i) + h * f(xi);
    % 校正步 (Corrector)
    y_imp_euler(i+1) = y_imp_euler(i) + (h/2) * (f(xi) + f(xi_next));
end

% 3. 打印结果表格
fprintf('Step n\t x_n\t Euler 数值解\t 改进 Euler 数值解\n');
fprintf('------------------------------------------------------------\n');
for i = 1:n
    fprintf('%d\t %.1f\t %.4f\t\t %.4f\n', i-1, x(i), y_euler(i), y_imp_euler(i));
end

% 4. 绘图对比
figure('Color', 'w');
plot(x, y_euler, 'ro-', 'DisplayName', '标准 Euler'); hold on;
plot(x, y_imp_euler, 'bs-', 'DisplayName', '改进 Euler');

% 为了对比，绘制 MATLAB 内置的高精度数值积分结果
integral_func = @(t) exp(t.^2);
x_fine = 0:0.01:1;
y_true = arrayfun(@(ax) integral(integral_func, 0, ax), x_fine);
plot(x_fine, y_true, 'k--', 'DisplayName', '参考精确曲线');

grid on;
xlabel('x');
ylabel('I(x)');
title(['步长 h = ', num2str(h), ' 时积分 I(x) = \int e^{t^2} dt 的计算']);
legend('Location', 'NorthWest');

% 输出最终结果
fprintf('\n在 x = 1 处的近似值：\n');
fprintf('标准 Euler 方法：%.4f\n', y_euler(end));
fprintf('改进 Euler 方法：%.4f\n', y_imp_euler(end));