function parallel_five_bar_leg_optimization()
%% 并联五连杆轮腿机器人跳跃高度PSO优化 - 增强版
% 参考专利: 《一种两轮足机器人跳跃姿态的力学分析及约束方法》(CN116738739A)
% 
% 优化策略:
%   1. 多次独立运行（5次），使用不同参数配置
%   2. 混沌初始化（Logistic映射）
%   3. 自适应PSO参数
%   4. 引入先验知识
%   5. 改进的收敛曲线可视化
%
% 跳跃阶段:
%   Stage1: 收腿，机体最低点，关节角度φ₁
%   Stage2: 伸腿，轮子离地前最后一刻，关节角度φ₂，机体位移h1
%   Stage3.1: 斜抛运动，保持φ₂，跳跃高度h2
%   Stage3.2: 空中收腿，关节角度φ₃，收腿高度h3
%
% 总跳跃高度: H = h2 + h3

clear; clc; close all;

%% ==================== 第一步：建立机器人运动学模型（符号推导）====================
fprintf('正在建立运动学模型...\n');
syms l1 l2 theta1 theta2 real

% 1.1 基本几何关系
% 以水平地面为x轴、以垂直于地的方向为y轴建立坐标系
% 关节A（原点）驱动 D 点
dy = l1 * cos(theta1);
dx = l1 * sin(theta1);

% 关节B（原点）驱动 E 点 (与A点重合，但控制不同的腿)
ey = l1 * cos(theta2);
ex = l1 * sin(theta2);

% 1.2 求解C点坐标（并联五连杆闭合方程）
% C点同时满足与D点距离为l2，与E点距离为l2
syms cx cy real
eq1 = (cx - dx)^2 + (cy - dy)^2 == l2^2;
eq2 = (cx - ex)^2 + (cy - ey)^2 == l2^2;
sol = solve([eq1, eq2], [cx, cy]);

% 选择物理可行的解（轮子在下方）
cx_sol = real(sol.cx(1));
cy_sol = real(sol.cy(1));

% 1.3 坐标变换（将坐标系旋转180度，使轮子朝下）
% 齐次变换矩阵
theta = 0;  % 不旋转，直接使用原坐标系
r = [cos(theta), -1*sin(theta), -1*cx_sol;
     sin(theta),  cos(theta),   -1*cy_sol;
     0,           0,            1];

% 对各点应用变换
d = [dx; dy; 1];
d_w = r * d;
dx_w = d_w(1);
dy_w = d_w(2);

e = [ex; ey; 1];
e_w = r * e;
ex_w = e_w(1);
ey_w = e_w(2);

c = [cx_sol; cy_sol; 1];
c_w = r * c;
cx_w = c_w(1);
cy_w = c_w(2);

a = [0; 0; 1];
a_w = r * a;
ax_w = a_w(1);
ay_w = a_w(2);

% 1.4 计算各杆质心坐标
% b: 小臂(小腿)质心位置参数, c_param: 大臂(大腿)质心位置参数
% MD=b*AD=b*l1, NC=c*CD=c*l2
syms b c_param real

% AD杆（大腿l1连杆）的质心M_D坐标
mx_d = (1-b)*dx_w + b * ax_w;
my_d = (1-b)*dy_w + b * ay_w;

% AE杆（大腿l1连杆）的质心M_E坐标
mx_e = (1-b)*ex_w + b * ax_w;
my_e = (1-b)*ey_w + b * ay_w;

% DC杆（小腿l2连杆）的质心N_D坐标
nx_d = (1-c_param)*cx_w + c_param * dx_w;
ny_d = (1-c_param)*cy_w + c_param * dy_w;

% EC杆（小腿l2连杆）的质心N_E坐标
nx_e = (1-c_param)*cx_w + c_param * ex_w;
ny_e = (1-c_param)*cy_w + c_param * ey_w;

fprintf('运动学模型建立完成。\n');

%% ==================== 第二步：推导跳跃高度函数（完整动力学）====================
fprintf('推导跳跃高度函数...\n');
syms phi1 phi2 phi3 t b_val c_val m_l1 m_l2 m_body m_wheel g real

% 2.1 定义三个关键状态的角度
% phi为关节电机与竖直向上方向的夹角
% 为保证机体只在竖直方向移动，规定2个关节电机输出的角度绝对值相等
% theta1 = phi, theta2 = -phi

% phi1时刻（Stage1: 缩腿最低点）
ay_phi1 = subs(ay_w, [theta1, theta2], [phi1, -1*phi1]);
dy_phi1 = subs(dy_w, [theta1, theta2], [phi1, -1*phi1]);
ey_phi1 = subs(ey_w, [theta1, theta2], [phi1, -1*phi1]);
my_d_phi1 = subs(my_d, [theta1, theta2], [phi1, -1*phi1]);
my_e_phi1 = subs(my_e, [theta1, theta2], [phi1, -1*phi1]);
ny_d_phi1 = subs(ny_d, [theta1, theta2], [phi1, -1*phi1]);
ny_e_phi1 = subs(ny_e, [theta1, theta2], [phi1, -1*phi1]);

% phi2时刻（Stage2: 离地瞬间）
ay_phi2 = subs(ay_w, [theta1, theta2], [phi2, -1*phi2]);
ax_phi2 = subs(ax_w, [theta1, theta2], [phi2, -1*phi2]);
dy_phi2 = subs(dy_w, [theta1, theta2], [phi2, -1*phi2]);
dx_phi2 = subs(dx_w, [theta1, theta2], [phi2, -1*phi2]);
ey_phi2 = subs(ey_w, [theta1, theta2], [phi2, -1*phi2]);
ex_phi2 = subs(ex_w, [theta1, theta2], [phi2, -1*phi2]);
my_d_phi2 = subs(my_d, [theta1, theta2], [phi2, -1*phi2]);
mx_d_phi2 = subs(mx_d, [theta1, theta2], [phi2, -1*phi2]);
my_e_phi2 = subs(my_e, [theta1, theta2], [phi2, -1*phi2]);
mx_e_phi2 = subs(mx_e, [theta1, theta2], [phi2, -1*phi2]);
ny_d_phi2 = subs(ny_d, [theta1, theta2], [phi2, -1*phi2]);
ny_e_phi2 = subs(ny_e, [theta1, theta2], [phi2, -1*phi2]);

% phi3时刻（Stage3.2: 最高点收腿）
ay_phi3 = subs(ay_w, [theta1, theta2], [phi3, -1*phi3]);
dy_phi3 = subs(dy_w, [theta1, theta2], [phi3, -1*phi3]);
ey_phi3 = subs(ey_w, [theta1, theta2], [phi3, -1*phi3]);

% 2.2 计算各阶段高度变化
% h1: A点位移量（Stage1到Stage2）
h1 = ay_phi2 - ay_phi1;

% h3: 主动收腿高度（Stage3.2）
h3 = ay_phi2 - ay_phi3;

% 2.3 各杆质心位移量（phi1到phi2）
% h_AD = M_Dy_phi2 - M_Dy_phi1 (AD杆质心位移)
h_AD = my_d_phi2 - my_d_phi1;
% h_AE = M_Ey_phi2 - M_Ey_phi1 (AE杆质心位移)  
h_AE = my_e_phi2 - my_e_phi1;
% h_DC = N_Dy_phi2 - N_Dy_phi1 (DC杆质心位移)
h_DC = ny_d_phi2 - ny_d_phi1;
% h_EC = N_Ey_phi2 - N_Ey_phi1 (EC杆质心位移)
h_EC = ny_e_phi2 - ny_e_phi1;

% 2.4 速度瞬心法求刚体各点速度
% Stage1到Stage2阶段，驱动轮不离开地面，C点固定
% 根据五连杆几何关系和两关节电机输出角度绝对值相等，机体沿竖直"滑轨"移动

% 速度瞬心图定义:
% CD杆为构件1; AD杆为构件2; AE杆为构件3; EC杆为构件4; 
% 机体为构件3; 固定支座和滑轨为构件6

% 瞬心P26的位置（AD杆与机体/滑轨的速度瞬心）
p26_y = ay_phi2;
p26_x = p26_y / (dy_phi2/dx_phi2);
AP26 = p26_x - ax_phi2;

% 计算各关键距离（用于速度瞬心法）
DP26 = real(sqrt((dx_phi2 - p26_x)^2 + (dy_phi2 - p26_y)^2));
MP26 = real(sqrt((mx_d_phi2 - p26_x)^2 + (my_d_phi2 - p26_y)^2));
AM_D = real(sqrt((ax_phi2 - mx_d_phi2)^2 + (ay_phi2 - my_d_phi2)^2));
AM_E = real(sqrt((ax_phi2 - mx_e_phi2)^2 + (ay_phi2 - my_e_phi2)^2));

% 瞬心P46的位置（AE杆与机体/滑轨的速度瞬心）
p46_y = ay_phi2;
p46_x = p46_y / (ey_phi2/ex_phi2);
AP46 = p46_x - ax_phi2;
EP46 = sqrt((ex_phi2 - p46_x)^2 + (ey_phi2 - p46_y)^2);
MP46 = sqrt((mx_e_phi2 - p46_x)^2 + (my_e_phi2 - p46_y)^2);

% 2.5 以机体线速度v_b为基准，构建各点速度关系
syms v_b real

% Vb × AP26 = ω (角速度)
% Vm = ω × AP26; V_AD = AP26/DP26 × Vb (式2.2.3-11)
w_AD = v_b / AP26;
v_wheel = 0;  % 轮子C点固定在地面
v_A = DP26 / AP26 * v_b;
v_N_dc = b_val * v_A;  % V_CD = c × V_m (式2.2.3-10)
v_M_ad = w_AD * MP26;

% 同理对AE、EC杆
% Vb × AP46 = ω; V_E = AP46/EP46 × Vb (式2.2.3-13)
w_AE = v_b / AP46;
v_E = EP46 / AP46 * v_b;
v_N_ec = b_val * v_E;  % V_CE = c × V_E (式2.2.3-14)
v_M_ae = w_AE * MP46;

% 2.6 计算各速度的竖直分量（用于动量和能量计算）
% 使用余弦定理求速度方向角
% cosM_P26_A = (AP26² + MP26² - AM_D²) / (2×AP26×MP26) (式2.2.3-21)
cosM_P26_A = (AP26^2 + MP26^2 - AM_D^2) / (2*AP26*MP26);
v_M_ad_y = v_M_ad * cosM_P26_A;

% V_AEy = V_AE × (AP46² + EP46² - AE²) / (AP46×EP46) (式2.2.3-23)
cosM_P46_A = (AP46^2 + MP46^2 - AM_E^2) / (2*AP46*MP46);
v_M_ae_y = v_M_ae * cosM_P46_A;

% V_CDy = V_CD × sin∠CDR = V_CD × Dx_phi2/l2 (式2.2.3-20)
v_N_dc_y = v_N_dc * dx_phi2 / l2;
% V_CEy = V_CE × Ex_phi2/l2 (式2.2.3-22)
v_N_ec_y = v_N_ec * ex_phi2 / l2;

% Vby = Vb (机体速度即为竖直方向速度)
v_b_y = v_b;
% Vwm = 0 (轮子固定)
v_w_y = 0;

% 2.7 机械能守恒分析（Stage1到Stage2）
% 关节电机以恒定扭矩τ做功: W = ∫τdφ = τ(φ₂ - φ₁) ... (式2.2.3-1)
% 两个关节电机总功: W = 2×t×(phi2-phi1)
W = 2*t*(phi2 - phi1);

% 应用机械能守恒定律（式2.2.3-2）:
% W = (m1×g×h_AD + m1×g×h_AE + m2×g×h_DC + m2×g×h_EC + mb×g×h1)
%   + (1/2×m1×v_AD² + 1/2×m1×v_AE² + 1/2×m2×v_DC² + 1/2×m2×v_EC² + 1/2×mb×vb²)

% 势能增加
PE = m_l1*g*h_AD + m_l1*g*h_AE + m_l2*g*h_DC + m_l2*g*h_EC + m_body*g*h1;

% 动能（在phi2时刻，只考虑竖直方向速度分量的平方）
KE = (1/2)*m_l1*v_M_ad^2 + (1/2)*m_l1*v_M_ae^2 + ...
     (1/2)*m_l2*v_N_dc^2 + (1/2)*m_l2*v_N_ec^2 + ...
     (1/2)*m_body*v_b^2;

% 能量方程: W = PE + KE
eq_energy = W == PE + KE;

% 求解机体速度v_b
v_b_solutions = solve(eq_energy, v_b);
% 选择正值解（向上运动）
v_b_func = v_b_solutions(2);  % 取正值

% 2.8 计算腾空高度h2（Stage2到Stage3.1）
% 在该阶段，机器人各部件相对静止，整体作斜抛运动
% 水平方向速度由驱动轮提供，竖直方向速度由Stage1到Stage2的关节电机输出扭矩提供

% 整体质量
mz = 2*m_l1 + 2*m_l2 + m_body + m_wheel;

% 根据动量守恒（式2.2.3-19）:
% mz×Vy = V_AEy×m1 + V_ADy×m1 + V_DCy×m2 + V_ECy×m2 + Vby×mb + Vwy×mw
syms v_y_total real
eq_momentum = mz*v_y_total == m_l1*v_M_ad_y + m_l1*v_M_ae_y + ...
                              m_l2*v_N_dc_y + m_l2*v_N_ec_y + ...
                              m_body*v_b_y + m_wheel*v_w_y;
% 代入v_b求解整体竖直速度
v_y_expr = solve(eq_momentum, v_y_total);
v_y_func = subs(v_y_expr, v_b, v_b_func);

% 利用机械能守恒（Stage2到Stage3.1）（式2.2.3-17, 18）:
% 1/2×mz×Vy² = mz×g×h2
% Vy = sqrt(2×g×h2)
% h2 = Vy²/(2×g)
h2 = v_y_func^2 / (2*g);

% 2.9 总跳跃高度
H = h2 + h3;

fprintf('跳跃高度函数推导完成。\n');

%% ==================== 第三步：建立约束条件 ====================
fprintf('建立约束条件...\n');

% 3.1 平均作用力约束（式2.2.4-1, 2.2.4-2）
% F_avg = W/h1 = τ(φ₂-φ₁)/h1
% 需要约束: 0 < F_avg ≤ (mz - mw)×g
F_avg = W / h1;

% 3.2 角速度约束
% Stage1到Stage2的时间t1（式2.2.4-4）:
% a1 = F_avg/mz
% t1 = (mz/F)×sqrt(2×g×h1)
% 平均角速度 ω̄₁ = (φ₂-φ₁)/t1 (式2.2.4-5)
a1 = F_avg / mz;
t1 = sqrt(2*g*h1) / a1;
w1_avg = (phi2 - phi1) / t1;

% Stage3.1到Stage3.2的时间t2:
% t2 = sqrt(2×h3/g) (自由落体时间的逆过程)
% 平均角速度 ω̄₂ = (φ₂-φ₃)/t2 (式2.2.4-6)
t2 = sqrt(2*h3/g);
w2_avg = (phi2 - phi3) / t2;

% 简化的作用力约束: F = 2×t/l1
F_simple = 2*t / l1;

fprintf('约束条件建立完成。\n');

%% ==================== 第四步：生成数值计算函数 ====================
fprintf('生成数值计算函数...\n');

% 4.1 固定参数赋值
m_l1_val = 0.5;        % 大腿质量 l1连杆 (kg)
m_l2_val = 0.6;        % 小腿质量 l2连杆 (kg)
m_body_val = 8.5;      % 机体质量 (kg)
m_wheel_val = 0.58;    % 轮子质量 (kg)
g_val = 9.81;          % 重力加速度 (m/s^2)
b_num = 0.5;           % 小腿质心位置参数
c_num = 0.85;          % 大腿质心位置参数

% 代入物理参数
H_func_pso = subs(H, [m_l1, m_l2, m_body, m_wheel, g, b_val, c_val, b, c_param], ...
                 [m_l1_val, m_l2_val, m_body_val, m_wheel_val, g_val, b_num, c_num, b_num, c_num]);
             
w1_pso = subs(w1_avg, [m_l1, m_l2, m_body, m_wheel, g, b_val, c_val, b, c_param], ...
                 [m_l1_val, m_l2_val, m_body_val, m_wheel_val, g_val, b_num, c_num, b_num, c_num]);
             
w2_pso = subs(w2_avg, [m_l1, m_l2, m_body, m_wheel, g, b_val, c_val, b, c_param], ...
                 [m_l1_val, m_l2_val, m_body_val, m_wheel_val, g_val, b_num, c_num, b_num, c_num]);
             
F_pso = subs(F_simple, [m_l1, m_l2, m_body, m_wheel, g, b_val, c_val, b, c_param], ...
                 [m_l1_val, m_l2_val, m_body_val, m_wheel_val, g_val, b_num, c_num, b_num, c_num]);

% 整体质量数值
mz_num = double(subs(mz, [m_l1, m_l2, m_body, m_wheel], ...
                    [m_l1_val, m_l2_val, m_body_val, m_wheel_val]));

% 4.2 创建数值计算函数句柄
H_func_handle = matlabFunction(H_func_pso, 'Vars', {l1, l2, phi1, phi2, phi3, t});
w1_handle = matlabFunction(w1_pso, 'Vars', {l1, l2, phi1, phi2, t});
w2_handle = matlabFunction(w2_pso, 'Vars', {l1, l2, phi1, phi2, phi3, t});
F_handle = matlabFunction(F_pso, 'Vars', {l1, l2, phi1, phi2, t});

fprintf('数值计算函数生成完成。\n');

%% ==================== 第五步：多运行增强PSO优化 ====================
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║           增强型PSO优化 - 多次独立运行策略                      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

% 5.1 优化参数设置
n_particles = 3000;      % 每次运行的粒子数量
max_iter = 800;          % 每次运行的最大迭代次数
dim = 6;                 % 优化变量维度 [l1, l2, phi1, phi2, phi3, t]
n_runs = 5;              % 独立运行次数

% 5.2 变量范围约束
lb = [0.1, 0.1, 60, 90, 60, 10];    % 下限 [l1, l2, phi1, phi2, phi3, t]
ub = [0.2, 0.3, 90, 150, 150, 15];  % 上限

% 5.3 约束阈值
w1_max = 8.4;          % 角速度约束1 (rad/s)
w2_max = 8.4;          % 角速度约束2 (rad/s)
F_max = (mz_num - m_wheel_val) * g_val;  % 作用力约束 (N)
l2_min_ratio = 1;      % 腿长比例约束 l2/l1 >= 1

% 5.4 先验知识（文档中的最优参数，可根据实际调整）
prior_knowledge = [0.15, 0.20, 75, 120, 90, 12.5];  % [l1, l2, phi1, phi2, phi3, t]

% 5.5 自适应参数配置（5组不同参数）
% 每次运行使用不同的PSO参数组合
pso_configs = struct();
pso_configs(1).w_range = [0.4, 0.9];   pso_configs(1).c1 = 1.49; pso_configs(1).c2 = 1.49; pso_configs(1).neighborhood_ratio = 0.20;
pso_configs(2).w_range = [0.3, 0.8];   pso_configs(2).c1 = 1.70; pso_configs(2).c2 = 1.30; pso_configs(2).neighborhood_ratio = 0.25;
pso_configs(3).w_range = [0.5, 1.0];   pso_configs(3).c1 = 1.30; pso_configs(3).c2 = 1.70; pso_configs(3).neighborhood_ratio = 0.30;
pso_configs(4).w_range = [0.2, 0.7];   pso_configs(4).c1 = 1.50; pso_configs(4).c2 = 1.60; pso_configs(4).neighborhood_ratio = 0.35;
pso_configs(5).w_range = [0.4, 0.95];  pso_configs(5).c1 = 1.60; pso_configs(5).c2 = 1.40; pso_configs(5).neighborhood_ratio = 0.40;

% 5.6 存储各运行结果
all_best_fitness = zeros(n_runs, 1);
all_best_params = zeros(n_runs, dim);
all_convergence_history = zeros(max_iter, n_runs);

fprintf('开始多次独立优化运行...\n');
fprintf('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

%% ==================== 第六步：执行多次独立运行 ====================
for run_idx = 1:n_runs
    fprintf('\n【运行 %d/%d】开始...\n', run_idx, n_runs);
    fprintf('参数配置: w∈[%.2f, %.2f], c1=%.2f, c2=%.2f, 邻域比例=%.2f\n', ...
        pso_configs(run_idx).w_range(1), pso_configs(run_idx).w_range(2), ...
        pso_configs(run_idx).c1, pso_configs(run_idx).c2, ...
        pso_configs(run_idx).neighborhood_ratio);
    
    % 获取当前运行的PSO参数
    w_min = pso_configs(run_idx).w_range(1);
    w_max = pso_configs(run_idx).w_range(2);
    c1 = pso_configs(run_idx).c1;
    c2 = pso_configs(run_idx).c2;
    neighborhood_ratio = pso_configs(run_idx).neighborhood_ratio;
    neighborhood_size = max(2, round(n_particles * neighborhood_ratio));
    
    % 6.1 混沌初始化（Logistic映射）
    fprintf('  使用Logistic混沌映射初始化粒子群...\n');
    particles = chaotic_initialization(n_particles, dim, lb, ub);
    
    % 第一次运行时，将先验知识作为初始粒子之一
    if run_idx == 1
        particles(1, :) = prior_knowledge;
        fprintf('  已添加先验知识粒子: [%.3f, %.3f, %.1f, %.1f, %.1f, %.2f]\n', prior_knowledge);
    end
    
    % 速度限制
    v_max = 0.2 * (ub - lb);
    v_min = -v_max;
    
    % 初始化速度（小随机值）
    velocities = 0.1 * (rand(n_particles, dim) - 0.5) .* (ub - lb);
    
    % 6.2 初始化最佳位置和适应度
    personal_best = particles;
    personal_best_fitness = -inf(n_particles, 1);
    global_best = particles(1, :);
    global_best_fitness = -inf;
    
    % 邻域最优（局部拓扑）
    local_best = particles;
    local_best_fitness = personal_best_fitness;
    
    % 历史最优记录
    stagnation_counter = 0;
    stagnation_threshold = 40;
    
    % 6.3 初始评估所有粒子
    for i = 1:n_particles
        [fitness, ~] = evaluate_particle(particles(i,:), H_func_handle, ...
            w1_handle, w2_handle, F_handle, w1_max, w2_max, l2_min_ratio, F_max);
        personal_best_fitness(i) = fitness;
        if fitness > global_best_fitness
            global_best_fitness = fitness;
            global_best = particles(i, :);
        end
    end
    local_best_fitness = personal_best_fitness;
    local_best = personal_best;
    
    fprintf('  初始最优高度: %.2f mm\n', global_best_fitness * 1000);
    
    % 6.4 记录历史
    fitness_history = zeros(max_iter, 1);
    
    % 6.5 主优化循环
    for iter = 1:max_iter
        % 计算自适应惯性权重（非线性递减）
        progress = iter / max_iter;
        w_pso = w_max - (w_max - w_min) * progress^1.5;  % 非线性递减
        
        iter_fitness = zeros(n_particles, 1);
        
        % 更新邻域最优（环形拓扑）
        for i = 1:n_particles
            % 获取邻域粒子索引
            neighbors = get_neighborhood_indices(i, n_particles, neighborhood_size);
            
            % 找出邻域中的最优粒子
            [~, best_neighbor_idx] = max(personal_best_fitness(neighbors));
            local_best(i, :) = personal_best(neighbors(best_neighbor_idx), :);
            local_best_fitness(i) = personal_best_fitness(neighbors(best_neighbor_idx));
        end
        
        % 评估并更新所有粒子
        for i = 1:n_particles
            p = particles(i, :);
            
            % 评估粒子
            [fitness, h] = evaluate_particle(p, H_func_handle, w1_handle, ...
                w2_handle, F_handle, w1_max, w2_max, l2_min_ratio, F_max);
            iter_fitness(i) = fitness;
            
            % 更新个体最佳
            if fitness > personal_best_fitness(i)
                personal_best_fitness(i) = fitness;
                personal_best(i, :) = p;
            end
            
            % 更新全局最佳
            if fitness > global_best_fitness
                global_best_fitness = fitness;
                global_best = p;
                stagnation_counter = 0;
            end
        end
        
        % 记录历史
        fitness_history(iter) = global_best_fitness;
        
        % 计算种群多样性
        diversity = mean(std(particles));
        
        % 检测停滞
        if iter > 1 && abs(fitness_history(iter) - fitness_history(iter-1)) < 1e-8
            stagnation_counter = stagnation_counter + 1;
        else
            stagnation_counter = 0;
        end
        
        % 更新粒子速度和位置
        for i = 1:n_particles
            r1 = rand(1, dim);
            r2 = rand(1, dim);
            r3 = rand(1, dim);
            
            % 使用邻域最优的PSO速度更新（LBEST拓扑）
            cognitive = c1 * r1 .* (personal_best(i, :) - particles(i, :));
            social_local = c2 * r2 .* (local_best(i, :) - particles(i, :));
            social_global = 0.5 * c2 * r3 .* (global_best - particles(i, :));  % 全局引导
            
            velocities(i, :) = w_pso * velocities(i, :) + cognitive + social_local + social_global;
            
            % 速度限制
            velocities(i, :) = max(velocities(i, :), v_min);
            velocities(i, :) = min(velocities(i, :), v_max);
            
            % 位置更新
            particles(i, :) = particles(i, :) + velocities(i, :);
            
            % 边界处理（反弹策略）
            for d = 1:dim
                if particles(i, d) < lb(d)
                    particles(i, d) = lb(d);
                    velocities(i, d) = -0.5 * velocities(i, d);
                elseif particles(i, d) > ub(d)
                    particles(i, d) = ub(d);
                    velocities(i, d) = -0.5 * velocities(i, d);
                end
            end
        end
        
        % 自适应变异（基于停滞和多样性）
        if stagnation_counter > stagnation_threshold || diversity < 0.01
            mutation_rate = 0.25;
            mutation_strength = 0.4;
            particles = adaptive_mutation(particles, lb, ub, mutation_rate, mutation_strength, global_best);
            stagnation_counter = 0;
        elseif mod(iter, 25) == 0
            mutation_rate = 0.08;
            mutation_strength = 0.1;
            particles = adaptive_mutation(particles, lb, ub, mutation_rate, mutation_strength, global_best);
        end
        
        % 混沌扰动（每60次迭代）
        if mod(iter, 60) == 0
            particles = chaotic_disturbance(particles, lb, ub, 0.12);
        end
        
        % 精英重置（每120次迭代）
        if mod(iter, 120) == 0
            [~, sorted_idx] = sort(iter_fitness);
            worst_n = round(0.08 * n_particles);
            for j = 1:worst_n
                idx = sorted_idx(j);
                particles(idx, :) = global_best + 0.25 * (rand(1, dim) - 0.5) .* (ub - lb);
                particles(idx, :) = max(min(particles(idx, :), ub), lb);
                velocities(idx, :) = 0.05 * (rand(1, dim) - 0.5) .* (ub - lb);
            end
        end
        
        % 定期输出进度
        if mod(iter, 200) == 0
            fprintf('  迭代 %d/%d: 当前最优 %.2f mm, 多样性 %.4f\n', ...
                iter, max_iter, global_best_fitness*1000, diversity);
        end
    end
    
    % 6.6 存储本次运行结果
    all_best_fitness(run_idx) = global_best_fitness;
    all_best_params(run_idx, :) = global_best;
    all_convergence_history(:, run_idx) = fitness_history;
    
    fprintf('  【运行 %d 完成】最优高度: %.4f m (%.2f mm)\n', run_idx, global_best_fitness, global_best_fitness*1000);
    fprintf('  最优参数: l1=%.4f, l2=%.4f, phi1=%.1f°, phi2=%.1f°, phi3=%.1f°, t=%.2f Nm\n', ...
        global_best(1), global_best(2), global_best(3), global_best(4), global_best(5), global_best(6));
end

%% ==================== 第七步：结果分析与显示 ====================
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║                      优化结果汇总                               ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

% 7.1 找出全局最优
[overall_best_fitness, best_run_idx] = max(all_best_fitness);
overall_best_params = all_best_params(best_run_idx, :);

% 7.2 计算统计信息
mean_fitness = mean(all_best_fitness);
std_fitness = std(all_best_fitness);
cv_fitness = std_fitness / mean_fitness * 100;  % 变异系数（%）

fprintf('\n━━━━━━━━━━ 各运行结果 ━━━━━━━━━━\n');
for i = 1:n_runs
    marker = '';
    if i == best_run_idx
        marker = ' ★ 最优';
    end
    fprintf('运行 %d: %.4f m (%.2f mm)%s\n', i, all_best_fitness(i), all_best_fitness(i)*1000, marker);
end

fprintf('\n━━━━━━━━━━ 统计分析 ━━━━━━━━━━\n');
fprintf('平均跳跃高度: %.4f m (%.2f mm)\n', mean_fitness, mean_fitness*1000);
fprintf('标准差:       %.6f m (%.4f mm)\n', std_fitness, std_fitness*1000);
fprintf('变异系数:     %.2f%%\n', cv_fitness);
fprintf('最优高度:     %.4f m (%.2f mm) [运行 %d]\n', overall_best_fitness, overall_best_fitness*1000, best_run_idx);
fprintf('最差高度:     %.4f m (%.2f mm)\n', min(all_best_fitness), min(all_best_fitness)*1000);

fprintf('\n━━━━━━━━━━ 全局最优参数 ━━━━━━━━━━\n');
fprintf('大腿长度 l1   = %.4f m\n', overall_best_params(1));
fprintf('小腿长度 l2   = %.4f m\n', overall_best_params(2));
fprintf('缩腿角度 phi1 = %.2f°\n', overall_best_params(3));
fprintf('离地角度 phi2 = %.2f°\n', overall_best_params(4));
fprintf('最高点角度 phi3 = %.2f°\n', overall_best_params(5));
fprintf('关节扭矩 t    = %.2f Nm\n', overall_best_params(6));
fprintf('腿长比例 l2/l1 = %.3f\n', overall_best_params(2)/overall_best_params(1));

% 7.3 优化稳定性评估
fprintf('\n━━━━━━━━━━ 优化稳定性评估 ━━━━━━━━━━\n');
if cv_fitness < 1
    stability_rating = '★★★★★ 非常稳定';
elseif cv_fitness < 3
    stability_rating = '★★★★☆ 较稳定';
elseif cv_fitness < 5
    stability_rating = '★★★☆☆ 一般稳定';
elseif cv_fitness < 10
    stability_rating = '★★☆☆☆ 稳定性较差';
else
    stability_rating = '★☆☆☆☆ 不稳定';
end
fprintf('稳定性评级: %s\n', stability_rating);

%% ==================== 第八步：可视化 ====================
% 8.1 创建Figure 1：所有运行的收敛曲线对比
figure('Name', '多运行收敛曲线对比', 'Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
colors = lines(n_runs);
hold on;
for i = 1:n_runs
    plot(1:max_iter, all_convergence_history(:, i) * 1000, ...
        'Color', colors(i, :), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('运行 %d', i));
end
xlabel('迭代次数', 'FontSize', 12);
ylabel('跳跃高度 (mm)', 'FontSize', 12);
title('所有运行的收敛曲线对比', 'FontSize', 14, 'FontWeight', 'bold');
legend('show', 'Location', 'southeast');
grid on;
ax = gca;
ax.YAxis.Exponent = 0;

% 8.2 最优运行的详细收敛曲线
subplot(1, 2, 2);
best_history = all_convergence_history(:, best_run_idx) * 1000;
plot(1:max_iter, best_history, 'b-', 'LineWidth', 2);
hold on;

% 标记关键点
[max_val, max_iter_idx] = max(best_history);
plot(max_iter_idx, max_val, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
text(max_iter_idx + 20, max_val, sprintf('最优: %.2f mm', max_val), ...
    'FontSize', 10, 'Color', 'r');

% 标记收敛阶段
convergence_threshold = 0.99 * max_val;
converged_iter = find(best_history >= convergence_threshold, 1);
if ~isempty(converged_iter)
    xline(converged_iter, '--g', sprintf('99%%收敛 @ iter %d', converged_iter), ...
        'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
end

xlabel('迭代次数', 'FontSize', 12);
ylabel('跳跃高度 (mm)', 'FontSize', 12);
title(sprintf('最优运行 (运行 %d) 收敛曲线', best_run_idx), 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ax = gca;
ax.YAxis.Exponent = 0;

% 8.3 创建Figure 2：各运行结果的柱状图对比
figure('Name', '各运行结果对比', 'Position', [150, 150, 800, 500]);

bar_data = all_best_fitness * 1000;
bar_colors = repmat([0.3, 0.6, 0.9], n_runs, 1);
bar_colors(best_run_idx, :) = [0.9, 0.3, 0.3];  % 最优运行用红色标记

b = bar(1:n_runs, bar_data, 'FaceColor', 'flat');
b.CData = bar_colors;

hold on;
% 添加均值线
yline(mean_fitness * 1000, '--k', sprintf('均值: %.2f mm', mean_fitness * 1000), ...
    'LineWidth', 2, 'LabelHorizontalAlignment', 'right');

% 添加数值标签
for i = 1:n_runs
    text(i, bar_data(i) + 0.5, sprintf('%.2f', bar_data(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end

xlabel('运行编号', 'FontSize', 12);
ylabel('最优跳跃高度 (mm)', 'FontSize', 12);
title('各运行最优结果对比', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% 添加图例说明
text(0.5, max(bar_data) * 0.95, sprintf('变异系数 CV = %.2f%%', cv_fitness), ...
    'FontSize', 11, 'Color', [0.5, 0.5, 0.5]);

% 8.4 创建Figure 3：参数分布箱线图
figure('Name', '参数分布分析', 'Position', [200, 200, 1200, 400]);

param_names = {'l_1 (m)', 'l_2 (m)', '\phi_1 (°)', '\phi_2 (°)', '\phi_3 (°)', '\tau (Nm)'};
for i = 1:dim
    subplot(1, dim, i);
    boxplot(all_best_params(:, i), 'Colors', 'b', 'Widths', 0.6);
    hold on;
    scatter(ones(n_runs, 1), all_best_params(:, i), 50, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
    ylabel(param_names{i}, 'FontSize', 11);
    title(sprintf('参数 %d', i), 'FontSize', 12);
    grid on;
end
sgtitle('各运行最优参数分布', 'FontSize', 14, 'FontWeight', 'bold');

fprintf('\n优化完成！请查看图形窗口获取详细可视化结果。\n');

end

%% ==================== 辅助函数：混沌初始化（Logistic映射）====================
function particles = chaotic_initialization(n_particles, dim, lb, ub)
    % 使用Logistic混沌映射生成初始种群
    % 公式: x_{n+1} = μ × x_n × (1 - x_n), μ = 4
    
    mu = 4;  % Logistic映射参数
    
    % 80%使用混沌映射，20%使用均匀分布增加覆盖性
    n_chaotic = round(0.8 * n_particles);
    n_uniform = n_particles - n_chaotic;
    
    % 混沌映射初始化
    chaotic_particles = zeros(n_chaotic, dim);
    
    for d = 1:dim
        % 为每个维度使用不同的初始值
        x = 0.1 + 0.1 * rand();  % 随机初始值，避免不动点
        
        % 丢弃前100个值以消除瞬态效应
        for k = 1:100
            x = mu * x * (1 - x);
        end
        
        % 生成混沌序列
        for i = 1:n_chaotic
            x = mu * x * (1 - x);
            % 映射到变量范围
            chaotic_particles(i, d) = lb(d) + x * (ub(d) - lb(d));
        end
    end
    
    % 均匀分布初始化（增加覆盖性）
    uniform_particles = lb + rand(n_uniform, dim) .* (ub - lb);
    
    % 合并
    particles = [chaotic_particles; uniform_particles];
    
    % 打乱顺序
    particles = particles(randperm(n_particles), :);
    
    % 确保在边界内
    particles = max(particles, lb);
    particles = min(particles, ub);
end

%% ==================== 辅助函数：获取邻域索引（环形拓扑）====================
function neighbors = get_neighborhood_indices(particle_idx, n_particles, neighborhood_size)
    % 获取指定粒子的邻域索引（环形拓扑）
    
    half_size = floor(neighborhood_size / 2);
    neighbors = zeros(neighborhood_size, 1);
    
    for j = 1:neighborhood_size
        offset = j - half_size - 1;
        neighbor_idx = mod(particle_idx - 1 + offset, n_particles) + 1;
        neighbors(j) = neighbor_idx;
    end
    
    neighbors = unique(neighbors);  % 去除重复
end

%% ==================== 辅助函数：粒子评估 ====================
function [fitness, h] = evaluate_particle(p, H_func_handle, w1_handle, ...
    w2_handle, F_handle, w1_max, w2_max, l2_min_ratio, F_max)
    % 提取参数
    l1_val = p(1); l2_val = p(2);
    phi1_val = deg2rad(p(3)); 
    phi2_val = deg2rad(p(4)); 
    phi3_val = deg2rad(p(5));
    t_val = p(6);
    
    % 计算目标函数值（跳跃高度）
    try
        h = H_func_handle(l1_val, l2_val, phi1_val, phi2_val, phi3_val, t_val);
        % 检查结果有效性
        if ~isreal(h) || isnan(h) || isinf(h) || h < 0 || h > 1
            h = 0;
        end
    catch
        h = 0;
    end
    
    % 计算约束条件
    try
        w1 = w1_handle(l1_val, l2_val, phi1_val, phi2_val, t_val);
        w2 = w2_handle(l1_val, l2_val, phi1_val, phi2_val, phi3_val, t_val);
        F_val = F_handle(l1_val, l2_val, phi1_val, phi2_val, t_val);
        l2_ratio = l2_val / l1_val;
    catch
        w1 = inf; w2 = inf; F_val = inf; l2_ratio = 0;
    end
    
    % 检查约束并计算惩罚项
    penalty = 0;
    
    % 角速度约束1
    if ~isreal(w1) || isnan(w1)
        penalty = penalty + 100;
    elseif abs(w1) > w1_max
        penalty = penalty + 100 * (abs(w1) - w1_max);
    end
    
    % 角速度约束2
    if ~isreal(w2) || isnan(w2)
        penalty = penalty + 100;
    elseif abs(w2) > w2_max
        penalty = penalty + 100 * (abs(w2) - w2_max);
    end
    
    % 作用力约束
    if ~isreal(F_val) || isnan(F_val) || isinf(F_val)
        penalty = penalty + 100;
    else
        if abs(F_val) > F_max
            penalty = penalty + 100 * (abs(F_val) - F_max);
        end
        if F_val < 0
            penalty = penalty + 100 * abs(F_val);
        end
    end
    
    % 腿长比例约束
    if l2_ratio < l2_min_ratio
        penalty = penalty + 100 * (l2_min_ratio - l2_ratio);
    end
    
    % 适应度计算（带惩罚项）
    fitness = h - penalty;
    
    % 确保输出合理
    if ~isreal(fitness) || isnan(fitness)
        fitness = -inf;
        h = 0;
    end
end

%% ==================== 辅助函数：混沌扰动 ====================
function particles = chaotic_disturbance(particles, lb, ub, strength)
    % 混沌映射（Logistic映射）扰动
    chaos = 0.7;
    mu = 4;  % Logistic参数
    
    for i = 1:size(particles, 1)
        if rand() < 0.25  % 25%概率扰动
            chaos = mu * chaos * (1 - chaos);
            perturbation = strength * (ub - lb) .* (chaos - 0.5);
            particles(i, :) = particles(i, :) + perturbation;
            particles(i, :) = max(min(particles(i, :), ub), lb);
        end
    end
end

%% ==================== 辅助函数：自适应变异 ====================
function particles = adaptive_mutation(particles, lb, ub, mutation_rate, mutation_strength, global_best)
    % 自适应变异操作
    
    n_particles = size(particles, 1);
    dim = size(particles, 2);
    
    for i = 1:n_particles
        if rand() < mutation_rate
            mutation_type = rand();
            
            if mutation_type < 0.3
                % 类型1：完全随机重置
                particles(i, :) = lb + rand(1, dim) .* (ub - lb);
                
            elseif mutation_type < 0.6
                % 类型2：高斯变异
                sigma = mutation_strength * (ub - lb);
                particles(i, :) = particles(i, :) + sigma .* randn(1, dim);
                
            elseif mutation_type < 0.8
                % 类型3：向全局最优方向变异
                direction = global_best - particles(i, :);
                step = mutation_strength * rand() * direction;
                particles(i, :) = particles(i, :) + step;
                
            else
                % 类型4：柯西变异
                cauchy_val = tan(pi * (rand(1, dim) - 0.5));
                cauchy_val = max(min(cauchy_val, 10), -10);
                particles(i, :) = particles(i, :) + mutation_strength * (ub - lb) .* cauchy_val * 0.1;
            end
            
            particles(i, :) = max(min(particles(i, :), ub), lb);
        end
    end
end
