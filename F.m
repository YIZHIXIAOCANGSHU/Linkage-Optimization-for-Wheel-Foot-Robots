%% 轮腿机器人 形变高度差Ay、形变速度v
clear;clc;close all;

%% 定义符号变量 集合参数
syms theta1 theta2 l1 l2 l3 k real

m = 15; % 單上质量(kg)
g = 9.81; % 重力加速度
r = 67.5; % 轮子半径(mm)

syms theta2 l1 l2 real
% 关节 H 点坐标
syms hx hy real
hx = l1 * cos(theta2);
hy = l1 * sin(theta2);

% 关节 K 点坐标，就是 J 点，用来求行程和行程速度(mm)
syms ky real
eq3 = (hx)^2 + (ky - hy)^2 == (l2)^2;
ky_sol = solve(eq3, ky);
for i = 1:length(ky_sol)
    current_sol = ky_sol(i);
    % 检查当前解是否为 a + b
    if simplify(current_sol - l1*sin(theta2)) == (l2 + l1*cos(theta2))^(1/2)*(l2 - l1*cos(theta2))^(1/2)
        rky_sol = current_sol * 1000 + r; % 单位换算成 mm 并加上轮子半径即可得到 h
        break; % 找到后退出循环
    end
end

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


%% 求最高点,最低点，行程
rky_max_num = double(subs(rky_sol, [l1, l2, l3, k, theta2], [l1_val, l2_val, l3_val, k_val, theta_max]));
rky_min_num = double(subs(rky_sol, [l1, l2, l3, k, theta2], [l1_val, l2_val, l3_val, k_val, theta_min]));
delta_y = rky_max_num - rky_min_num;
fprintf('l1=%.f,l2=%.f,l3=%.f 时，最高点=%.2f,最低点=%.2f，行程=%.2f\n', l1_mm, l2_mm, l3_mm, rky_max_num, rky_min_num, delta_y);

theta_range = min_ang:0.1:max_ang;
theta_range_val = deg2rad(theta_range);
k_rky = subs(rky_sol, [l1, l2, l3, k], [l1_val, l2_val, l3_val, k_val]);
k_rky_num = matlabFunction(k_rky); % 将符号表达式转换为数值函数
f_vals = k_rky_num(theta_range_val); % 计算 f(x) 的值

% 线性拟合
p = polyfit(theta_range, f_vals, 1);
slope = p(1);
intercept = p(2);
fit_line = polyval(p, theta_range);

fprintf('形变的平均速率为: %f\n', slope);

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