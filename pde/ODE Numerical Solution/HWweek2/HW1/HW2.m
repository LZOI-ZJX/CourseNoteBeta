%% 数值方法对比：改进 Euler 与 RK4
clear; clc;

% 1. 参数设置
h = 0.5;                % 步长
x_end = 2;              % 范围
x = 0:h:x_end;
n = length(x);

% 初始化结果向量
y_euler = zeros(1, n);  % 标准 Euler (一阶)
y_imp_euler = zeros(1, n); % 改进 Euler (二阶)
y_rk4 = zeros(1, n);    % RK4 (四阶)

% 初值条件
y_euler(1) = 0;
y_imp_euler(1) = 0;
y_rk4(1) = 0;

% 定义导数函数 f(x, y)
f = @(x, y) 1 - (2*x*y) / (1 + x^2);
% 定义精确解函数
y_exact_func = @(x) (x.*(3 + x.^2)) ./ (3.*(1 + x.^2));

% 2. 迭代计算
for i = 1:(n-1)
    xi = x(i);
    yi_e = y_euler(i);
    yi_ie = y_imp_euler(i);
    yi_rk = y_rk4(i);
    
    % --- 标准 Euler 方法 ---
    y_euler(i+1) = yi_e + h * f(xi, yi_e);
    
    % --- 改进 Euler 方法 (Predictor-Corrector) ---
    yp = yi_ie + h * f(xi, yi_ie); % 预估步
    y_imp_euler(i+1) = yi_ie + (h/2) * (f(xi, yi_ie) + f(xi+h, yp)); % 校正步
    
    % --- 四阶 RK 方法 ---
    k1 = f(xi, yi_rk);
    k2 = f(xi + h/2, yi_rk + h*k1/2);
    k3 = f(xi + h/2, yi_rk + h*k2/2);
    k4 = f(xi + h, yi_rk + h*k3);
    y_rk4(i+1) = yi_rk + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

% 3. 计算精确值
y_exact = y_exact_func(x);

% 4. 打印结果对比表
fprintf('x_n\t\t精确解\t\tEuler\t\t改进Euler\tRK4\n');
fprintf('----------------------------------------------------------------------\n');
for i = 1:n
    fprintf('%.1f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\n', ...
        x(i), y_exact(i), y_euler(i), y_imp_euler(i), y_rk4(i));
end

% 5. 绘图对比
figure('Color', 'w');
x_fine = 0:0.01:x_end;
plot(x_fine, y_exact_func(x_fine), 'k-', 'LineWidth', 1.5, 'DisplayName', 'Exact Solution'); hold on;
plot(x, y_euler, 'ro--', 'DisplayName', 'Standard Euler (1st)');
plot(x, y_imp_euler, 'gs--', 'DisplayName', 'Improved Euler (2nd)');
plot(x, y_rk4, 'bd--', 'DisplayName', 'RK4 (4th)');

grid on; xlabel('x'); ylabel('y');
title('不同数值方法精度对比 (h=0.5)');
legend('Location', 'SouthEast');