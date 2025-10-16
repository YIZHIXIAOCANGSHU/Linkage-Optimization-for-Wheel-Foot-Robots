%% 并联五连杆轮腿机器人关节电机扭矩负载计算
clear; clc; close all;
%% --- 参数定义 ---
% 机器人基本参数
m = 15;        % 机器人质量 (kg)
g = 9.81;      % 重力加速度 (m/s^2)

% 连杆长度参数（单位：米）
l1_val = 0.145; % 主动连杆长度
l2_val = 0.27;  % 从动连杆长度
l3_val = 0.1;   % 基座宽度

% 关节角度参数（以度给出，内部转换为弧度）
a_ang = 0;      % A 点角度（度）
b_ang = 50;     % B 点角度（度）


%% --- 符号变量与几何关系 ---
% 定义关节角度和连杆长度为实数符号变量
syms theta1 theta2 l1 l2 l3 real

% A 点驱动 D 点坐标计算
% D点是第一个连杆的末端，基于角度theta1和长度l1计算坐标
Dx = l1 * cos(theta1) + l3/2;
Dy = l1 * sin(theta1);

% B 点（位于 (l3,0)）驱动 E 点坐标计算
% E点是第二个连杆的末端，基于角度theta2和长度l1计算坐标
Ex = l1 * cos(theta2) - l3/2;
Ey = l1 * sin(theta2);

% C点满足与D、E两点距离均为l2的约束条件
% 建立两个圆的方程，求解它们的交点
syms Cx Cy real

eq1 = (Cx - Dx)^2 + (Cy - Dy)^2 == l2^2;  % C点到D点距离为l2
eq2 = (Cx - Ex)^2 + (Cy - Ey)^2 == l2^2;  % C点到E点距离为l2

% 求解两个圆的交点（通常会有两个解）
sol = solve([eq1, eq2], [Cx, Cy]);

%% --- 参数预处理 ---
% 角度参数转换为弧度
theta1_val = deg2rad(180 - a_ang); % 转换为标准坐标系角度
theta2_val = deg2rad(b_ang);

%% --- 选择物理解 ---
% 将符号解代入数值参数进行计算
substitution_vars = [l1 l2 l3 theta1 theta2];
substitution_vals = [l1_val l2_val l3_val theta1_val theta2_val];

% 计算两个解的坐标值
cx1_num = double(subs(sol.Cx(1), substitution_vars, substitution_vals));
cy1_num = double(subs(sol.Cy(1), substitution_vars, substitution_vals));
cx2_num = double(subs(sol.Cx(2), substitution_vars, substitution_vals));
cy2_num = double(subs(sol.Cy(2), substitution_vars, substitution_vals));

% 检查解的有效性
if cy1_num <= 0 && cy2_num <= 0
    error('两个解的y值都小于等于0，不符合物理意义！')
end

% 选择 y 值较大的有效解（通常为上方交点）
% 如果两个解的y值都大于0，则选择y值较大的
% 如果只有一个解的y值大于0，则选择该解
if cy1_num > 0 && cy2_num > 0
    % 两个解都有效，选择y值较大的
    if cy1_num >= cy2_num
        idx = 1;
    else
        idx = 2;
    end
elseif cy1_num > 0
    % 只有第一个解有效
    idx = 1;
else
    % 只有第二个解有效
    idx = 2;
end

% 取出被选择的符号解（作为雅可比的基准）
Cx_sym = simplify(sol.Cx(idx));
Cy_sym = simplify(sol.Cy(idx));

%% --- 雅可比矩阵（对 theta1, theta2 求偏导）---
% 构建雅可比矩阵 J = [dCx/dθ1  dCx/dθ2; dCy/dθ1  dCy/dθ2]
% 雅可比矩阵描述了末端点C的速度与关节角度变化率的关系
J_sym = jacobian([Cx_sym; Cy_sym], [theta1, theta2]);
% 在代入几何参数和角度后得到数值雅可比矩阵
J_num = double(subs(J_sym, substitution_vars, substitution_vals));

%% --- 关节扭矩计算 ---
% 末端受力（示例：假设地面反作用力为 m*g/2 向上，横向为0）
% 根据力平衡原理，计算关节需要产生的扭矩
F = [0; m * g / 2]; % 单位 N，2x1 向量 (Fx; Fy)

% 将末端力映射到关节扭矩：tau = J^T * F
% 使用雅可比矩阵的转置将末端力转换为关节扭矩
tau = J_num.' * F;      % 单位 N·m
tau_mm = tau * 1000;    % 单位 N·mm

%% --- 输出结果 ---
fprintf('A 点关节角为 %d 度，B 点关节角为 %d 度\n', a_ang, b_ang);
disp('=============== 并联五连杆雅可比矩阵 (数值) ===============');
disp(J_num);
disp('=============== 并联五连杆关节扭矩 ===============');
fprintf('τ1 = %.6f N·m (%.3f N·mm)\n', tau(1), tau_mm(1));
fprintf('τ2 = %.6f N·m (%.3f N·mm)\n', tau(2), tau_mm(2));