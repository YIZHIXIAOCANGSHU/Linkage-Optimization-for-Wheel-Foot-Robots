%% 偏置并联五连杆轮腿机器人k值变化与关键轴系负载
clear; clc; close all;

%% ========== 参数定义 ==========
%% 物理参数
m = 15;      % 质量(kg)
g = 9.81;    % 重力加速度
l1_val = 0.1;   % 连杆长度(m)
l2_val = 0.15;  % 连杆长度(m)
l4_val = 0.04;  % 新增距离参数(m)

%% 关节角度参数
a_ang = 30;  
b_ang = 30;  
theta1_val = deg2rad(a_ang);
theta2_val = deg2rad(180 - b_ang);

%% 符号变量定义
syms theta1 theta2 l1 l2 l3 l4 k real
syms Cx Cy real
syms Hx Hy real
syms Jx Jy real
syms Fx Fy real
syms Ex Ey real
syms Dx Dy real

%% ========== 遍历k值计算 ==========
k_values = 0.01:0.01:1;  % k值范围从0.01到1，步进0.01 (避免k=0时的除零错误)
F1_values = zeros(size(k_values));
F2_values = zeros(size(k_values));

for i = 1:length(k_values)
    k_val = k_values(i);

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
    substitution_vars = [l1 l2 l4 k theta1 theta2];
    substitution_vals = [l1_val l2_val l4_val k_val theta1_val theta2_val];
    
    % 计算两个解的坐标值
    Cx1_num = double(subs(sol.Cx(1), substitution_vars, substitution_vals));
    Cy1_num = double(subs(sol.Cy(1), substitution_vars, substitution_vals));
    Cx2_num = double(subs(sol.Cx(2), substitution_vars, substitution_vals));
    Cy2_num = double(subs(sol.Cy(2), substitution_vars, substitution_vals));
    
    % 检查解的有效性
    if Cy1_num <= 0 && Cy2_num <= 0
        warning('k=%.2f时，两个解的y值都小于等于0，不符合物理意义！', k_val);
        % 将无效值设为NaN
        F1_values(i) = NaN;
        F2_values(i) = NaN;
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
    
    %% d_JH J点到H点的水平距离 
    d_JH = abs(Hx - Jx); 
    
    %% d_CE 计算E点到DC直线的距离
    % 先计算DC向量
    DC_x = Dx - Cx_sol;
    DC_y = Dy - Cy_sol;
    
    % 计算DE向量
    DE_x = Dx - Ex;
    DE_y = Dy - Ey;
    
    % 利用向量叉积计算点到直线距离：d = |DC × DE| / |DC|
    cross_product = abs(DC_x * DE_y - DC_y * DE_x);
    DC_magnitude = sqrt(DC_x^2 + DC_y^2);
    d_CE = cross_product / DC_magnitude;
    
    %% d_GH
    % 计算F点坐标（在CE方向上，从E点延申l4距离）
    % 向量CE
    CE_x = Ex - Cx_sol;
    CE_y = Ey - Cy_sol;
    % CE向量的单位向量
    CE_magnitude = sqrt(CE_x^2 + CE_y^2);
    CE_unit_x = CE_x / CE_magnitude;
    CE_unit_y = CE_y / CE_magnitude;

    % F点坐标（从E点沿CE方向延伸l4）
    Fx = Ex + l4 * CE_unit_x;
    Fy = Ey + l4 * CE_unit_y;
    
    % EH向量
    EH_x = Hx - Ex;
    EH_y = Hy - Ey;
    
    % FG直线方程（过F点且平行于EH）
    % 直线方向向量与EH相同
    FG_direction_x = EH_x;
    FG_direction_y = EH_y;
    % FG直线方程可以表示为通过点F且方向向量为(EH_x, EH_y)的直线
    % Ax + By + C = 0 形式
    % A = FG_direction_y, B = -FG_direction_x, C = FG_direction_x * Fy - FG_direction_y * Fx
    A = FG_direction_y;
    B = -FG_direction_x;
    C = FG_direction_x * Fy - FG_direction_y * Fx;
    
    
    % 点到直线距离公式: d = |Ax + By + C| / sqrt(A^2 + B^2)
    d_GH = abs(A * Hx + B * Hy + C) / sqrt(A^2 + B^2);
    
    F=m * g / 2;

    % 计算F1和F2
    F1 = (d_JH/d_GH) * F;
    F2 = (d_JH/d_CE) * F;
    
    % 计算F1和F2的数值解
    % 先检查各个分母是否为零，避免除零错误
    
    % 检查F1 = (d_JH/d_GH) * F 中的分母 d_GH 是否为零
    % d_GH = abs(A * Hx + B * Hy + C) / sqrt(A^2 + B^2)
    % 需要先检查 sqrt(A^2 + B^2) 是否为零
    denominator_dGH = sqrt(A^2 + B^2);
    denominator_dGH_val = double(subs(denominator_dGH, substitution_vars, substitution_vals));
    if abs(denominator_dGH_val) < eps || isinf(denominator_dGH_val) || isnan(denominator_dGH_val)
        F1_values(i) = NaN;
    else
        % 计算d_GH的值
        d_GH_val = double(subs(d_GH, substitution_vars, substitution_vals));
        if abs(d_GH_val) < eps || isinf(d_GH_val) || isnan(d_GH_val)
            F1_values(i) = NaN;
        else
            F1_values(i) = double(subs(F1, substitution_vars, substitution_vals));
        end
    end
    
    % 检查F2 = (d_JH/d_CE) * F 中的分母 d_CE 是否为零
    % d_CE = cross_product / DC_magnitude
    % 需要先检查 DC_magnitude 是否为零
    DC_magnitude_val = double(subs(DC_magnitude, substitution_vars, substitution_vals));
    if abs(DC_magnitude_val) < eps || isinf(DC_magnitude_val) || isnan(DC_magnitude_val)
        F2_values(i) = NaN;
    else
        % 计算d_CE的值
        d_CE_val = double(subs(d_CE, substitution_vars, substitution_vals));
        if abs(d_CE_val) < eps || isinf(d_CE_val) || isnan(d_CE_val)
            F2_values(i) = NaN;
        else
            % CE_magnitude 在计算CE_unit时也会用到，检查是否为零
            CE_magnitude_val = double(subs(CE_magnitude, substitution_vars, substitution_vals));
            if abs(CE_magnitude_val) < eps || isinf(CE_magnitude_val) || isnan(CE_magnitude_val)
                % CE_magnitude为零会影响Fx,Fy的计算，进而可能影响d_GH的计算，这里设F2为NaN
                F2_values(i) = NaN;
            else
                % CE_unit_x = CE_x / CE_magnitude 和 CE_unit_y = CE_y / CE_magnitude 中检查CE_magnitude
                F2_values(i) = double(subs(F2, substitution_vars, substitution_vals));
            end
        end
    end

end

%% ========== 可视化 ==========
figure;
plot(k_values, F1_values, 'b-', 'LineWidth', 1.5);
hold on;
plot(k_values, F2_values, 'r--', 'LineWidth', 1.5);
xlabel('k值');
ylabel('负载力 (N)');
title('F1和F2随k值的变化曲线');
legend('F1', 'F2');
grid on;
xlim([min(k_values), max(k_values)]);
ylim([0, max(max(F1_values(~isnan(F1_values))), max(F2_values(~isnan(F2_values)))) * 1.1]);

% 显示一些统计信息
fprintf('F1最大值: %.2f N\n', max(F1_values(~isnan(F1_values))));
fprintf('F1最小值: %.2f N\n', min(F1_values(~isnan(F1_values))));
fprintf('F2最大值: %.2f N\n', max(F2_values(~isnan(F2_values))));
fprintf('F2最小值: %.2f N\n', min(F2_values(~isnan(F2_values))));