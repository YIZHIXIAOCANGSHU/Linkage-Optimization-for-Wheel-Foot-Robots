%%  串联腿机器人的扭矩分析

clear; clc; close all;

%% ========== 参数定义 ==========
%% 物理参数
m = 15;      % 质量(kg)
g = 9.81;    % 重力加速度
l1_val = 0.1;   % 连杆长度(m)
l2_val = 0.15;  % 连杆长度(m)


%% 关节角度参数
a_ang = 30;  
b_ang = 30;  
theta1_val = deg2rad(180 - a_ang);
theta2_val = deg2rad(b_ang);

%% 符号变量定义
syms theta1 theta2 l1 l2 l3 k real
syms Cx Cy real

%% ========== 正向运动学计算 ==========

Dx =  l1 * cos(theta1);
Dy =  l1 * sin(theta1);

%% D点坐标数值计算
Dx_num = double(subs(Dx, [l1, theta1], [l1_val, theta1_val]));
Dy_num = double(subs(Dy, [l1, theta1], [l1_val, theta1_val]));

fprintf('D点坐标: (%.4f, %.4f)\n', Dx_num, Dy_num);

%% 关于C点坐标，解方程组，符号解
eq1 = (Cx - Dx)^2 + (Cy - Dy)^2 == (l2)^2;
eq2 = (Dy -Cy)/(Dx - Cx) == tan(theta2);
sol = solve([eq1, eq2], [Cx, Cy]);

%% 提取解
solCx = sol.Cx;
solCy = sol.Cy;

%% C点坐标数值计算
Cx_solutions = double(subs(solCx, [l1, l2, theta1, theta2], [l1_val, l2_val, theta1_val, theta2_val]));
Cy_solutions = double(subs(solCy, [l1, l2, theta1, theta2], [l1_val, l2_val, theta1_val, theta2_val]));

fprintf('C点可能的坐标:\n');
for i = 1:length(Cx_solutions)
    fprintf('  解%d: (%.4f, %.4f)\n', i, Cx_solutions(i), Cy_solutions(i));
end

%% 选择合理的C点坐标（在D点左侧的解，即x坐标小于Dx_num的解）
% 查找x坐标小于D点x坐标的解
valid_indices = find(Cx_solutions < Dx_num);
if isempty(valid_indices)
    error('没有找到在D点左侧的C点解');
end

% 选择第一个满足条件的解
idx = valid_indices(1);

Cx_num = Cx_solutions(idx);
Cy_num = Cy_solutions(idx);

%% 存储对应的选择好的符号解
Cx_sol_selected = solCx(idx);
Cy_sol_selected = solCy(idx);

fprintf('选择的C点坐标: (%.4f, %.4f)\n', Cx_num, Cy_num);

%% ========== 雅可比矩阵与力矩计算 ==========
%% 对theta1 和 theta2 求偏导，使用选定的解
J = [diff(Cx_sol_selected, theta1), diff(Cx_sol_selected, theta2);
     diff(Cy_sol_selected, theta1), diff(Cy_sol_selected, theta2)];

%% 雅可比矩阵计算
J_num = subs(J, [l1, l2], [l1_val, l2_val]);
J_final = double(subs(J_num, [theta1, theta2], [theta1_val, theta2_val]));

%% 显示雅可比矩阵维度信息
fprintf('雅可比矩阵维度: %d x %d\n', size(J_final, 1), size(J_final, 2));

%% 关节力矩计算
F = [0; -m * g/2]; % 作用在C点的重力
tau_final = J_final' * F; % 关节力矩(N·m)
tau_final_mm = tau_final * 1000; % 转换为 N·mm

fprintf('A点关节力为%.f度，B点关节力为%.f度\n', a_ang, b_ang);
disp('==============================');
disp('J = ');
disp(J_final);
disp('==============================');
fprintf('t1 = %.2f N·m (%.2f N·mm)\n', tau_final(1), tau_final_mm(1));
fprintf('t2 = %.2f N·m (%.2f N·mm)\n', tau_final(2), tau_final_mm(2));