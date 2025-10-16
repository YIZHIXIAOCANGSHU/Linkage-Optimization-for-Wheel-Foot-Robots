%% 偏置并联五连杆轮腿机器人关节电机扭矩负载计算
clear; clc; close all;

%% ========== 参数定义 ==========
%% 物理参数
m = 15;      % 质量(kg)
g = 9.81;    % 重力加速度
l1_val = 0.1;   % 连杆长度(m)
l2_val = 0.15;  % 连杆长度(m)
k_val = 0.5;    % 系数 k < 1

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
sol = solve([eq1, eq2], [Cx, Cy]);

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
    error('两个解的y值都小于等于0，不符合物理意义！')
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
fprintf('J点坐标为 (%.2f, %.5f)\n', x_num, y_num);

%% 解算C点坐标
% 直接使用之前选择的解，而不是重新求解
cx_num = double(subs(Cx_sol, [l1, l2, k, theta1, theta2], [l1_val, l2_val, k_val, theta1_val, theta2_val]));
cy_num = double(subs(Cy_sol, [l1, l2, k, theta1, theta2], [l1_val, l2_val, k_val, theta1_val, theta2_val]));

fprintf('C点坐标为 (%.5f, %.5f)\n', cx_num, cy_num);

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

fprintf('A点关节力为%.f度，B点关节力为%.f度\n', a_ang, b_ang);
disp('==============================');
disp('J = ');
disp(J_final);
disp('==============================');
fprintf('t1 = %.2f N·m (%.2f N·mm)\n', tau_final(1), tau_final_mm(1));
fprintf('t2 = %.2f N·m (%.2f N·mm)\n', tau_final(2), tau_final_mm(2));