% 1D Heat Equation: Comparison of Finite Difference Schemes with Plotting
clear; clc; close all; % 添加了 close all 清除历史图像

% --- 基本参数设置 ---
h = 0.1;               % 空间步长
x = 0:h:1;             % 空间网格节点
N = length(x) - 1;     % 空间区间等分数 (N=10)
T_end = 0.1;           % 计算终止时间

lambdas = [0.1, 0.5];  % 两种网格比情况

% 循环处理两种不同的 \lambda
for L_idx = 1:length(lambdas)
    lambda = lambdas(L_idx);
    tau = lambda * h^2;         % 时间步长
    M = round(T_end / tau);     % 时间步数
    
    fprintf('==============================================================\n');
    fprintf('计算情况: lambda = %.1f, h = %.1f, tau = %.4f, 总步数 = %d\n', lambda, h, tau, M);
    fprintf('==============================================================\n');
    
    % --- 精确解 ---
    u_exact = exp(-pi^2 * T_end) * sin(pi * x);
    
    % --- 初始化初值 ---
    u0 = sin(pi * x); % t=0 时的初值
    
    %% 1. 向前差分格式 (Forward Euler)
    u_FE = u0;
    for k = 1:M
        u_next = zeros(1, N+1);
        for i = 2:N
            u_next(i) = u_FE(i) + lambda * (u_FE(i+1) - 2*u_FE(i) + u_FE(i-1));
        end
        u_FE = u_next;
    end
    
    %% 2. 向后差分格式 (Backward Euler)
    main_diag = (1 + 2*lambda) * ones(N-1, 1);
    off_diag = -lambda * ones(N-2, 1);
    A_BE = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);
    
    u_BE = u0;
    for k = 1:M
        b = u_BE(2:N)';
        u_inner = A_BE \ b;
        u_BE(2:N) = u_inner';
    end
    
    %% 3. Crank-Nicolson 格式
    main_diag_A = (1 + lambda) * ones(N-1, 1);
    off_diag_A = -lambda/2 * ones(N-2, 1);
    A_CN = diag(main_diag_A) + diag(off_diag_A, 1) + diag(off_diag_A, -1);
    
    main_diag_B = (1 - lambda) * ones(N-1, 1);
    off_diag_B = lambda/2 * ones(N-2, 1);
    B_CN = diag(main_diag_B) + diag(off_diag_B, 1) + diag(off_diag_B, -1);
    
    u_CN = u0;
    for k = 1:M
        b = B_CN * u_CN(2:N)';
        u_inner = A_CN \ b;
        u_CN(2:N) = u_inner';
    end
    
    %% 4. Du Fort-Frankel 格式
    u_prev = u0; 
    u_curr = exp(-pi^2 * tau) * sin(pi * x); % t = tau 时的精确解
    
    if M == 1
        u_DF = u_curr;
    else
        for k = 2:M
            u_next = zeros(1, N+1);
            for i = 2:N
                u_next(i) = ((1 - 2*lambda)*u_prev(i) + 2*lambda*(u_curr(i+1) + u_curr(i-1))) / (1 + 2*lambda);
            end
            u_prev = u_curr;
            u_curr = u_next;
        end
        u_DF = u_curr;
    end
    
    %% --- 输出比较结果 ---
    fprintf('%-4s | %-10s | %-10s | %-10s | %-10s | %-10s\n', 'x', 'Exact', 'Forward', 'Backward', 'Crank-N', 'Du Fort-F');
    fprintf('--------------------------------------------------------------\n');
    for i = 2:N
        fprintf('%.1f  | %10.6f | %10.6f | %10.6f | %10.6f | %10.6f\n', ...
            x(i), u_exact(i), u_FE(i), u_BE(i), u_CN(i), u_DF(i));
    end
    
    fprintf('\n【绝对误差最大值】\n');
    fprintf('向前差分   : %e\n', max(abs(u_exact - u_FE)));
    fprintf('向后差分   : %e\n', max(abs(u_exact - u_BE)));
    fprintf('Crank-Nic  : %e\n', max(abs(u_exact - u_CN)));
    fprintf('Du Fort-F  : %e\n\n', max(abs(u_exact - u_DF)));
    
    %% --- 绘图模块 ---
    figure('Name', sprintf('结果对比 (Lambda = %.1f)', lambda), 'Position', [100+L_idx*50, 100+L_idx*50, 700, 500]);
    
    % 绘制精确解 (黑色实线)
    plot(x, u_exact, 'k-', 'LineWidth', 2); hold on;
    
    % 绘制各数值解 (采用不同颜色和标记)
    plot(x, u_FE, 'r--o', 'LineWidth', 1.5, 'MarkerSize', 6);
    plot(x, u_BE, 'b-.*', 'LineWidth', 1.5, 'MarkerSize', 6);
    plot(x, u_CN, 'g:s',  'LineWidth', 1.5, 'MarkerSize', 6);
    plot(x, u_DF, 'm-.d', 'LineWidth', 1.5, 'MarkerSize', 6);
    
    % 图形修饰
    xlabel('空间位置 x', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('温度 u(x,t)', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('数值格式对比 (t = 0.1, \\lambda = %.1f)', lambda), 'FontSize', 14);
    legend('精确解', '向前差分 (Explicit)', '向后差分 (Implicit)', 'Crank-Nicolson', 'Du Fort-Frankel', 'Location', 'best');
    grid on;
    
    % 设置坐标轴范围让边缘更美观
    xlim([0 1]);
    
    hold off;
end