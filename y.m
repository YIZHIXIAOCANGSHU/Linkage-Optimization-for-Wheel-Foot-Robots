%% 形变高度差Ay、形变速度v表达式
clear;clc;close all;
clearvars;

%% 定义符号变量 集合参数
syms theta1 theta2 l1 l2 k real

m = 15; % 簧上质量(kg)
g = 9.81; % 重力加速度

%% 连杆长度
l1_mm = 145; % mm
l2_mm = 270; % mm
l3_mm = 0; % mm
l1_val = l1_mm / 1000; % m
l2_val = l2_mm / 1000; % m
l3_val = l3_mm / 1000; % m
k_val = 0.5; % 系数, k<1

%% 设置关节转角范围 (单位: °)
min_ang = 0;     % 最小角度，例如 0°
max_ang = 60;    % 最大角度，例如 60°

% 转换为弧度制
theta_min = deg2rad(min_ang);
theta_max = deg2rad(max_ang);

% 关节 H 点坐标
syms Hx Hy real
Hx = l1 * cos(theta2);
Hy = l1 * sin(theta2);

% 关节 K 点坐标，就是 J 点，用来求行程和行程速度(mm)
syms Kx Ky real
Kx_sol = 0;
eq3 = (Kx-Hx)^2 + (Ky - Hy)^2 == (l2)^2;

% 在一个循环中同时计算所有角度的解并选择有效解
theta_range = min_ang:0.1:max_ang;
theta_range_val = deg2rad(theta_range);
f_vals = zeros(size(theta_range_val)); % 预分配数组

% 初始化最大值和最小值
rky_max_num = -inf;
rky_min_num = inf;

for i = 1:length(theta_range_val)
    theta_curr = theta_range_val(i);
    
    % 求解Ky，临时关闭警告信息
    warning('off', 'symbolic:solve:UnsupportedConditions');
    warning('off', 'symbolic:solve:WarnIfParams');
    Ky_sol = solve(eq3, Ky, 'IgnoreAnalyticConstraints', true);
    warning('on', 'symbolic:solve:UnsupportedConditions');
    warning('on', 'symbolic:solve:WarnIfParams');

    
    % 计算两个解的坐标值（使用cell数组进行符号替换）
    Ky1_num = double(subs(Ky_sol(1), {l1, l2, Kx, theta2}, {l1_val, l2_val, Kx_sol, theta_curr}));
    Ky2_num = double(subs(Ky_sol(2), {l1, l2, Kx, theta2}, {l1_val, l2_val, Kx_sol, theta_curr}));
    
    % 检查解的有效性
    if Ky1_num <= 0 && Ky2_num <= 0
        error('两个解的y值都小于等于0，不符合物理意义！')
    end
    
    % 选择 y 值较大的有效解（通常为上方交点）
    % 如果两个解的y值都大于0，则选择y值较大的
    % 如果只有一个解的y值大于0，则选择该解
    if Ky1_num > 0 && Ky2_num > 0
        % 两个解都有效，选择y值较大的
        if Ky1_num >= Ky2_num
            idx = 1;
        else
            idx = 2;
        end
    elseif Ky1_num > 0
        % 只有第一个解有效
        idx = 1;
    else
        % 只有第二个解有效
        idx = 2;
    end
    
    % 取出被选择的符号解并计算数值
    Ky_selected = Ky_sol(idx);

    rky_curr = double(subs(Ky_selected * 1000, {l1, l2, Kx, theta2}, {l1_val, l2_val, Kx_sol, theta_curr}));
    
    f_vals(i) = rky_curr;
    
    % 更新最大值和最小值
    if rky_curr > rky_max_num
        rky_max_num = rky_curr;
    end
    if rky_curr < rky_min_num
        rky_min_num = rky_curr;
    end
end

delta_y = rky_max_num - rky_min_num;

% 线性拟合
p = polyfit(theta_range, f_vals, 1);
slope = p(1);
intercept = p(2);
fit_line = polyval(p, theta_range);

% 对当前输出终端进行清空
clc;

% 输出形变的平均速率
fprintf('形变的平均速率为: %f\n', slope);
% 输出形变的最大值和最小值以及差值
fprintf('形变的最大值为: %f\n', rky_max_num);
fprintf('形变的最小值为: %f\n', rky_min_num);
fprintf('形变的差值为: %f\n', delta_y);

%% 可视化结果
figure;
% 绘制原始数据
plot(theta_range, f_vals, 'b-', 'LineWidth', 1.5);
hold on;
plot(theta_range, fit_line, 'r--', 'LineWidth', 1.5);
xlabel('关节角度 \theta (^{\circ})');
ylabel('行程 (mm)');
title('行程 h vs 关节角度');
legend('实际行程', '线性拟合', 'Location', 'best');
grid on;
axis tight;

% 添加拟合方程文本
eqn_text = sprintf('拟合直线: y = %.4f x + %.4f', slope, intercept);
text(mean(theta_range), mean(f_vals) - 20, eqn_text, ...
    'FontSize', 10, 'BackgroundColor', 'white');