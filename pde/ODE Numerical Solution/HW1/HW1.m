%% Euler 方法求解初值问题
clear; clc;

% 1. 参数设置
h = 0.5;                % 步长
x_end = 2;              % 范围终点
x = 0:h:x_end;          % 离散点 x_n
n = length(x);
y_euler = zeros(1, n);  % 初始化数值解向量
y_euler(1) = 0;         % 初值 y(0) = 0

% 定义导数函数 f(x, y)
f = @(x, y) 1 - (2*x*y) / (1 + x^2);

% 定义精确解函数 y(x)
y_exact_func = @(x) (x.*(3 + x.^2)) ./ (3.*(1 + x.^2));

% 2. Euler 迭代计算 (显式方法)
for i = 1:(n-1)
    y_euler(i+1) = y_euler(i) + h * f(x(i), y_euler(i));
end

% 3. 计算精确解
y_exact = y_exact_func(x);

% 4. 打印结果表格
fprintf('x_n\t\t数值解 (Euler)\t精确解\t\t绝对误差\n');
fprintf('------------------------------------------------------------\n');
for i = 1:n
    err = abs(y_exact(i) - y_euler(i));
    fprintf('%.1f\t\t%.5f\t\t%.5f\t\t%.5f\n', x(i), y_euler(i), y_exact(i), err);
end

% 5. 绘图对比
figure('Color', 'w');
x_fine = 0:0.01:x_end; % 用于绘制平滑的精确解曲线
plot(x_fine, y_exact_func(x_fine), 'b-', 'LineWidth', 1.5); hold on;
plot(x, y_euler, 'ro--', 'MarkerSize', 8, 'LineWidth', 1);

grid on;
xlabel('x');
ylabel('y');
title('Euler方法数值解与精确解对比');
legend('精确解', ['Euler 方法 (h=', num2str(h), ')'], 'Location', 'SouthEast');