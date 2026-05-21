%% 四阶 RK 方法求解三阶 ODE
clear; clc;

h = 0.2;
x = 0:h:3;
N = length(x);
Y = zeros(3, N);
Y(:,1) = [1; 2; 0]; % 初值 [y; y'; y'']

% 定义方程组 F(x, Y)
F = @(x, Y) [Y(2); ...
             Y(3); ...
             exp(x) - 2*Y(3) + Y(2) + 2*Y(1)];

% 迭代计算
for n = 1:N-1
    xn = x(n);
    Yn = Y(:,n);
    
    k1 = F(xn, Yn);
    k2 = F(xn + h/2, Yn + h/2*k1);
    k3 = F(xn + h/2, Yn + h/2*k2);
    k4 = F(xn + h, Yn + h*k3);
    
    Y(:,n+1) = Yn + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
end

% 结果提取
y_val = Y(1,:);

% 绘图
plot(x, y_val, 'r-o', 'LineWidth', 1.5);
grid on; xlabel('x'); ylabel('y');
title('四阶 RK 方法求解结果 (h=0.2)');