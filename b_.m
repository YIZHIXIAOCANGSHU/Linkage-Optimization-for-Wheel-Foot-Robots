%% 偏置并联五连杆轮腿机器人k值变化与关节电机扭矩负载
clear; clc; close all;

%% ========== 参数定义 ==========
%% 物理参数
m = 15;      % 质量(kg)
g = 9.81;    % 重力加速度
l1_val = 0.1;   % 连杆长度(m)
l2_val = 0.15;  % 连杆长度(m)
% k_val = 0.5;    % 系数 k < 1  移除此行，因为k将作为遍历变量

%% 关节角度参数
a_ang = 30;  
b_ang = 30;  
theta1_val = deg2rad(a_ang);
theta2_val = deg2rad(180 - b_ang);

%% 符号变量定义
syms theta1 theta2 l1 l2 l3 k real
syms Cx Cy real
syms Hx Hy real
syms Jx Jy real

%% ========== 遍历k值计算 ==========
k_values = 0.01:0.01:1;  % k值范围从0.01到1，步进0.01 (避免k=0时的除零错误)
t1_values = zeros(size(k_values));
t2_values = zeros(size(k_values));

for i = 1:length(k_values)
    k_val = k_values(i);
    
    %% ========== 正向运动学计算 ==========
    %% 节点A（原点）驱动E点
    Ex = k * l1 * cos(theta1);
    Ey = k * l1 * sin(theta1);
    %% 节点A（原点）驱动D点
    % k < 1
    Dx = k * l1 * cos(theta2);
    Dy = k * l1 * sin(theta2);
    
    %% 关于C点坐标，解方程组，符号解
    eq1 = (Cx - Dx)^2 + (Cy - Dy)^2 == (k * l2)^2;
    eq2 = (Cx - Ex)^2 + (Cy - Ey)^2 == (k * l2)^2;
    
    % 使用ReturnConditions参数来获取解及其条件，这样可以避免警告
    % 这种方法比忽略解析约束更严格，能获得完整的解信息同时避免警告
    [sol] = solve([eq1, eq2], [Cx, Cy]);
    
    %% 选择物理可行的解
    % 将符号解代入数值参数进行计算
    substitution_vars = [l1 l2 k theta1 theta2];
    substitution_vals = [l1_val l2_val k_val theta1_val theta2_val];
    
    % 计算两个解的坐标值
    Cx1_num = double(subs(sol.Cx(1), substitution_vars, substitution_vals));
    Cy1_num = double(subs(sol.Cy(1), substitution_vars, substitution_vals));
    Cx2_num = double(subs(sol.Cx(2), substitution_vars, substitution_vals));
    Cy2_num = double(subs(sol.Cy(2), substitution_vars, substitution_vals));
    
    % 检查解的有效性
    if Cy1_num <= 0 && Cy2_num <= 0
        warning('k=%.2f时，两个解的y值都小于等于0，不符合物理意义！', k_val);
        continue;
    end
    
    % 选择 y 值较大的有效解（通常为上方交点）
    % 如果两个解的y值都大于0，则选择y值较大的
    % 如果只有一个解的y值大于0，则选择该解
    if Cy1_num > 0 && Cy2_num > 0
        % 两个解都有效，选择y值较大的
        if Cy1_num >= Cy2_num
            idx = 1;
        else
            idx = 2;
        end
    elseif Cy1_num > 0
        % 只有第一个解有效
        idx = 1;
    else
        % 只有第二个解有效
        idx = 2;
    end
    
    % 取出被选择的符号解
    Cx_sol = simplify(sol.Cx(idx));
    Cy_sol = simplify(sol.Cy(idx));
    
    %% 关于H点坐标
    Hx = l1 * cos(theta1);
    Hy = l1 * sin(theta1);
    
    %% 关于J点坐标
    Jx = Hx + 1/k * (Cx_sol - Ex);
    Jy = Hy + 1/k * (Cy_sol - Ey);
    
    %% 计算J点坐标
    x_num = double(subs(Jx, [l1, l2, k, theta1, theta2], [l1_val, l2_val, k_val, theta1_val, theta2_val]));
    y_num = double(subs(Jy, [l1, l2, k, theta1, theta2], [l1_val, l2_val, k_val, theta1_val, theta2_val]));
    % fprintf('J点坐标为 (%.2f, %.5f)\n', x_num, y_num);
    
    %% 解算C点坐标
    % 直接使用之前选择的解，而不是重新求解
    cx_num = double(subs(Cx_sol, [l1, l2, k, theta1, theta2], [l1_val, l2_val, k_val, theta1_val, theta2_val]));
    cy_num = double(subs(Cy_sol, [l1, l2, k, theta1, theta2], [l1_val, l2_val, k_val, theta1_val, theta2_val]));
    
    % fprintf('C点坐标为 (%.5f, %.5f)\n', cx_num, cy_num);
    
    %% ========== 雅可比矩阵与力矩计算 ==========
    %% 对theta1 和 theta2 求偏导
    J = [diff(Jx, theta1), diff(Jx, theta2);
         diff(Jy, theta1), diff(Jy, theta2)];
    
    %% 雅可比矩阵计算
    J_num = subs(J, [l1, l2, k], [l1_val, l2_val, k_val]);
    J_final = double(subs(J_num, [theta1, theta2], [theta1_val, theta2_val]));
    
    %% 关节力矩计算
    F = [0; -m * g / 2]; % 每连杆上作用力
    tau_final = J_final' * F; % 关节力矩(N·m)
    tau_final_mm = tau_final * 1000; % 转换为 N·mm
    
    % 保存结果
    t1_values(i) = tau_final(1);
    t2_values(i) = tau_final(2);
    
    % 显示部分结果
    if mod(i, 20) == 1 || i == length(k_values)  % 每20个点显示一次进度
        fprintf('k=%.2f: t1=%.2f N·m, t2=%.2f N·m\n', k_val, tau_final(1), tau_final(2));
    end
end

%% ========== 结果可视化 ==========
figure;
plot(k_values, t1_values, 'b-', 'LineWidth', 1.5); hold on;
plot(k_values, t2_values, 'r-', 'LineWidth', 1.5);
xlabel('k 值');
ylabel('关节扭矩 (N·m)');
title('关节扭矩随k值变化');
legend('t1 (A点关节扭矩)', 't2 (B点关节扭矩)', 'Location', 'best');
grid on;

%% ========== 结果分析 ==========
fprintf('\n========== 结果分析 ==========\n');
fprintf('k值范围: %.2f 到 %.2f\n', min(k_values), max(k_values));
fprintf('t1 (A点关节扭矩):\n');
fprintf('  最大值: %.2f N·m (k=%.2f时)\n', max(t1_values), k_values(find(t1_values==max(t1_values), 1)));
fprintf('  最小值: %.2f N·m (k=%.2f时)\n', min(t1_values), k_values(find(t1_values==min(t1_values), 1)));
fprintf('  平均值: %.2f N·m\n', mean(t1_values));
fprintf('t2 (B点关节扭矩):\n');
fprintf('  最大值: %.2f N·m (k=%.2f时)\n', max(t2_values), k_values(find(t2_values==max(t2_values), 1)));
fprintf('  最小值: %.2f N·m (k=%.2f时)\n', min(t2_values), k_values(find(t2_values==min(t2_values), 1)));
fprintf('  平均值: %.2f N·m\n', mean(t2_values));

% 检查t1和t2是否随k变化
if abs(max(t1_values) - min(t1_values)) < 1e-10
    fprintf('\n结论: t1 (A点关节扭矩) 在k值变化时基本保持不变\n');
else
    fprintf('\n结论: t1 (A点关节扭矩) 随k值变化而变化\n');
end

if abs(max(t2_values) - min(t2_values)) < 1e-10
    fprintf('结论: t2 (B点关节扭矩) 在k值变化时基本保持不变\n');
else
    fprintf('结论: t2 (B点关节扭矩) 随k值变化而变化\n');
end