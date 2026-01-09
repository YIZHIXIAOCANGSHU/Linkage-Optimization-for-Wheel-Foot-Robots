%% 两轮足机器人跳跃高度单点验证计算
% 基于专利 CN 116738739 A 的数学模型
% 使用表2中的最优参数进行单点计算验证
% 目标: 验证能否得到 280.5mm 的跳跃总高度
%
% 重要说明：
% - 跳跃总高度 H = H1(斜抛高度) + H2(空中主动缩腿高度)
% - 表2中的"地面伸腿高度"是起跳过程中重心上升的几何高度，不是H1
% - H1 = H - H2 = 280.5 - 259.9 = 20.6mm

clear; clc; close all;

fprintf('============================================\n');
fprintf('  两轮足机器人跳跃高度 - 单点验证计算\n');
fprintf('============================================\n\n');

%% ========== 步骤1: 机器人参数初始化 (表1) ==========
fprintf('>>> 步骤1: 加载机器人基础参数 (表1)\n');
fprintf('--------------------------------------------\n');

% 质量参数 (kg)
m1 = 1.7 * 2;       % 轮子质量 (两个轮子)
m2 = 0.8 * 2;       % 小腿质量 (两条小腿)
m3 = 1.6 * 2;       % 大腿质量 (两条大腿)
m4 = 2.8;           % 躯体质量
mz = 11;            % 整体质量

% 几何参数 (m)
lz = 0.4;           % 大小腿总长度 400mm
r1 = 0.0445;        % 踝关节电机半径 44.5mm
r2 = 0.046;         % 髋关节电机半径 46mm

% 质心位置比例
b = 0.5;            % 小腿质心比例
c = 0.8;            % 大腿质心比例

% 关节参数
M2 = 33.5;          % 膝关节输出力矩 (Nm)
w_max = 21;         % 膝关节最大转速 (rad/s)

% 物理常数
g = 9.8;            % 重力加速度 (m/s^2)

fprintf('轮子质量 m1 = %.1f kg (1.7×2)\n', m1);
fprintf('小腿质量 m2 = %.1f kg (0.8×2)\n', m2);
fprintf('大腿质量 m3 = %.1f kg (1.6×2)\n', m3);
fprintf('躯体质量 m4 = %.1f kg\n', m4);
fprintf('整体质量 mz = %.0f kg\n', mz);
fprintf('大小腿总长度 lz = %.0f mm\n', lz * 1000);
fprintf('膝关节输出力矩 M2 = %.1f Nm\n', M2);
fprintf('膝关节最大转速 ωmax = %.0f rad/s\n', w_max);
fprintf('\n');

%% ========== 步骤2: 设定最优参数 (表2) ==========
fprintf('>>> 步骤2: 加载最优参数 (表2)\n');
fprintf('--------------------------------------------\n');

% 小腿长度占比
k = 0.4675;

% 膝关节角度 (度 -> 弧度)
eps1_deg = 147.3934;    % 地面缩腿膝关节角度
eps2_deg = 69.4081;     % 刚离地膝关节角度
eps3_deg = 147.3934;    % 飞行阶段膝关节角度

% 髋关节角度 (度 -> 弧度)
delta1_deg = 84.1274;   % 地面缩腿髋关节角度
delta2_deg = 36.8412;   % 刚离地髋关节角度
delta3_deg = 114.8207;  % 飞行阶段髋关节角度

% 转换为弧度
eps1 = deg2rad(eps1_deg);
eps2 = deg2rad(eps2_deg);
eps3 = deg2rad(eps3_deg);
delta1 = deg2rad(delta1_deg);
delta2 = deg2rad(delta2_deg);
delta3 = deg2rad(delta3_deg);

% 计算腿长
l1 = k * lz;            % 小腿长度
l2 = (1 - k) * lz;      % 大腿长度

fprintf('小腿长度占比 k = %.4f\n', k);
fprintf('小腿长度 l1 = %.2f mm\n', l1 * 1000);
fprintf('大腿长度 l2 = %.2f mm\n', l2 * 1000);
fprintf('\n');
fprintf('地面缩腿膝关节角度 ε1 = %.4f°\n', eps1_deg);
fprintf('刚离地膝关节角度 ε2 = %.4f°\n', eps2_deg);
fprintf('飞行阶段膝关节角度 ε3 = %.4f°\n', eps3_deg);
fprintf('地面缩腿髋关节角度 δ1 = %.4f°\n', delta1_deg);
fprintf('刚离地髋关节角度 δ2 = %.4f°\n', delta2_deg);
fprintf('飞行阶段髋关节角度 δ3 = %.4f°\n', delta3_deg);
fprintf('\n');

%% ========== 步骤3: 计算主动缩腿高度 H2 (式5) ==========
fprintf('>>> 步骤3: 计算主动缩腿高度 H2 (式5)\n');
fprintf('--------------------------------------------\n');

% 式(5): H2 = l1[cos(ε2-δ2) - cos(ε3-δ3)] + l2(cosδ2 - cosδ3)
term_H2_1 = l1 * (cos(eps2 - delta2) - cos(eps3 - delta3));
term_H2_2 = l2 * (cos(delta2) - cos(delta3));
H2 = term_H2_1 + term_H2_2;

fprintf('H2 = l1×[cos(ε2-δ2) - cos(ε3-δ3)] + l2×(cosδ2 - cosδ3)\n');
fprintf('   = %.4f × [cos(%.4f) - cos(%.4f)] + %.4f × (cos(%.4f) - cos(%.4f))\n', ...
    l1*1000, rad2deg(eps2-delta2), rad2deg(eps3-delta3), l2*1000, delta2_deg, delta3_deg);
fprintf('第一项: l1部分 = %.4f mm\n', term_H2_1 * 1000);
fprintf('第二项: l2部分 = %.4f mm\n', term_H2_2 * 1000);
fprintf('>>> 主动缩腿高度 H2 = %.4f mm\n', H2 * 1000);
fprintf('\n');

%% ========== 步骤4: 计算"地面伸腿高度" ==========
fprintf('>>> 步骤4: 计算"地面伸腿高度"（起跳过程重心上升高度）\n');
fprintf('--------------------------------------------\n');

% 地面伸腿阶段，机器人从(ε1,δ1)变化到(ε2,δ2)
% 重心上升高度 = l1[cos(ε2-δ2) - cos(ε1-δ1)] + l2[cos(δ2) - cos(δ1)]
h_ground_extend = l1 * (cos(eps2 - delta2) - cos(eps1 - delta1)) + ...
                  l2 * (cos(delta2) - cos(delta1));

fprintf('地面伸腿高度 = l1×[cos(ε2-δ2) - cos(ε1-δ1)] + l2×[cos(δ2) - cos(δ1)]\n');
fprintf('cos(ε1-δ1) = cos(%.4f°) = %.6f\n', rad2deg(eps1-delta1), cos(eps1-delta1));
fprintf('cos(ε2-δ2) = cos(%.4f°) = %.6f\n', rad2deg(eps2-delta2), cos(eps2-delta2));
fprintf('cos(δ1) = cos(%.4f°) = %.6f\n', delta1_deg, cos(delta1));
fprintf('cos(δ2) = cos(%.4f°) = %.6f\n', delta2_deg, cos(delta2));
fprintf('>>> 地面伸腿高度 = %.4f mm (文档值: 222.1 mm)\n', h_ground_extend * 1000);
fprintf('\n');

%% ========== 步骤5: 根据约束条件计算平均作用力 F̄ ==========
fprintf('>>> 步骤5: 根据约束条件反推平均作用力 F̄\n');
fprintf('--------------------------------------------\n');

% 约束条件: 0 < F̄ ≤ (mz - m1)×g
% 文档表2给出的结果暗示 F̄ 取约束上限附近的值
F_upper = (mz - m1) * g;
fprintf('平均作用力约束上限: F̄_max = (mz-m1)×g = (%.1f-%.1f)×%.1f = %.2f N\n', ...
    mz, m1, g, F_upper);

% 从式(22)计算理论F̄值（用于参考）
sin_diff1 = sin(eps1 - delta1);
sin_diff2 = sin(eps2 - delta2);
F_bar_theory = 2 * M2 / (l1 * (sin_diff2 - sin_diff1)) * log(sin_diff2 / sin_diff1);
fprintf('式(22)理论计算值: F̄_theory = %.2f N\n', F_bar_theory);
fprintf('注意: 理论值超出约束，实际运行时会被电机限制\n');

% 使用约束上限作为实际F̄值（这是优化结果的隐含条件）
F_bar = F_upper;
fprintf('>>> 实际使用的平均作用力 F̄ = %.2f N (约束上限)\n', F_bar);
fprintf('\n');

%% ========== 步骤6: 计算G系数 (式18) ==========
fprintf('>>> 步骤6: 计算G系数 (式18)\n');
fprintf('--------------------------------------------\n');

% 式(18): G = (m2*b+m3)*cos²δ2 + (1/2)*(m3*(c-1)-m2*b)*cot(ε2)*sin(2δ2) + m3*c*sin²δ2 + m4
cos_d2 = cos(delta2);
sin_d2 = sin(delta2);
cot_e2 = cot(eps2);
sin_2d2 = sin(2 * delta2);

term_G1 = (m2*b + m3) * cos_d2^2;
term_G2 = 0.5 * (m3*(c-1) - m2*b) * cot_e2 * sin_2d2;
term_G3 = m3 * c * sin_d2^2;
term_G4 = m4;

G = term_G1 + term_G2 + term_G3 + term_G4;

fprintf('G = (m2×b+m3)×cos²δ2 + 0.5×(m3×(c-1)-m2×b)×cot(ε2)×sin(2δ2) + m3×c×sin²δ2 + m4\n');
fprintf('cosδ2 = %.6f, sinδ2 = %.6f\n', cos_d2, sin_d2);
fprintf('cot(ε2) = %.6f, sin(2δ2) = %.6f\n', cot_e2, sin_2d2);
fprintf('第一项: (m2×b+m3)×cos²δ2 = %.6f\n', term_G1);
fprintf('第二项: 0.5×(m3×(c-1)-m2×b)×cot(ε2)×sin(2δ2) = %.6f\n', term_G2);
fprintf('第三项: m3×c×sin²δ2 = %.6f\n', term_G3);
fprintf('第四项: m4 = %.6f\n', term_G4);
fprintf('>>> G = %.6f\n', G);
fprintf('\n');

%% ========== 步骤7: 计算K系数 (式25) ==========
fprintf('>>> 步骤7: 计算K系数 (式25)\n');
fprintf('--------------------------------------------\n');

% 式(25): K = m2*b² × (cos²δ2/sin²ε2) + m3*cos²δ2 × [1 + (c*tanδ2 - (1-c)*cotε2)²] + m4
sin_e2 = sin(eps2);
tan_d2 = tan(delta2);

term_K1 = m2 * b^2 * (cos_d2^2 / sin_e2^2);
inner_term = c * tan_d2 - (1 - c) * cot_e2;
term_K2 = m3 * cos_d2^2 * (1 + inner_term^2);
term_K3 = m4;

K = term_K1 + term_K2 + term_K3;

fprintf('K = m2×b²×(cos²δ2/sin²ε2) + m3×cos²δ2×[1+(c×tanδ2-(1-c)×cotε2)²] + m4\n');
fprintf('sinε2 = %.6f, tanδ2 = %.6f\n', sin_e2, tan_d2);
fprintf('内部项: c×tanδ2-(1-c)×cotε2 = %.6f\n', inner_term);
fprintf('第一项: m2×b²×(cos²δ2/sin²ε2) = %.6f\n', term_K1);
fprintf('第二项: m3×cos²δ2×[1+(...)²] = %.6f\n', term_K2);
fprintf('第三项: m4 = %.6f\n', term_K3);
fprintf('>>> K = %.6f\n', K);
fprintf('\n');

%% ========== 步骤8: 计算T和Q (式25) ==========
fprintf('>>> 步骤8: 计算T和Q (式25)\n');
fprintf('--------------------------------------------\n');

% 式(25):
% T = cosδ2 - cosδ1
% Q = cos(ε2-δ2) - cos(ε1-δ1)
T = cos(delta2) - cos(delta1);
Q = cos(eps2 - delta2) - cos(eps1 - delta1);

fprintf('T = cosδ2 - cosδ1 = cos(%.4f°) - cos(%.4f°)\n', delta2_deg, delta1_deg);
fprintf('  = %.6f - %.6f = %.6f\n', cos(delta2), cos(delta1), T);
fprintf('>>> T = %.6f\n', T);
fprintf('\n');
fprintf('Q = cos(ε2-δ2) - cos(ε1-δ1) = cos(%.4f°) - cos(%.4f°)\n', ...
    rad2deg(eps2-delta2), rad2deg(eps1-delta1));
fprintf('  = %.6f - %.6f = %.6f\n', cos(eps2-delta2), cos(eps1-delta1), Q);
fprintf('>>> Q = %.6f\n', Q);
fprintf('\n');

%% ========== 步骤9: 计算斜抛运动高度 H1 (式24, 17) ==========
fprintf('>>> 步骤9: 计算斜抛运动高度 H1\n');
fprintf('--------------------------------------------\n');

% 使用式(24)计算v4²，然后用式(17)反推H1
% 式(24): v4² = (2/K) × {[F̄-(m2b+m3+m4)g]l1Q + [F̄-(m3c+m4)g]l2T}
% 式(17): v4 = mz × sqrt(2gH1) / G  =>  H1 = (v4 × G / mz)² / (2g)

bracket1 = F_bar - (m2*b + m3 + m4) * g;
bracket2 = F_bar - (m3*c + m4) * g;

fprintf('F̄ = %.2f N (使用约束上限)\n', F_bar);
fprintf('(m2×b+m3+m4)×g = %.2f N\n', (m2*b+m3+m4)*g);
fprintf('(m3×c+m4)×g = %.2f N\n', (m3*c+m4)*g);
fprintf('\n');
fprintf('第一个括号: F̄-(m2×b+m3+m4)×g = %.2f - %.2f = %.4f N\n', ...
    F_bar, (m2*b+m3+m4)*g, bracket1);
fprintf('第二个括号: F̄-(m3×c+m4)×g = %.2f - %.2f = %.4f N\n', ...
    F_bar, (m3*c+m4)*g, bracket2);

% 式(24)计算v4²
term1 = bracket1 * l1 * Q;
term2 = bracket2 * l2 * T;
v4_squared = (2/K) * (term1 + term2);
v4 = sqrt(v4_squared);

fprintf('\n');
fprintf('式(24): v4² = (2/K) × {[...]×l1×Q + [...]×l2×T}\n');
fprintf('第一项: [...]×l1×Q = %.4f × %.4f × %.6f = %.6f\n', bracket1, l1, Q, term1);
fprintf('第二项: [...]×l2×T = %.4f × %.4f × %.6f = %.6f\n', bracket2, l2, T, term2);
fprintf('v4² = (2/%.4f) × (%.6f + %.6f) = %.6f m²/s²\n', K, term1, term2, v4_squared);
fprintf('>>> 离地速度 v4 = %.4f m/s\n', v4);

% 式(17)计算H1
% v4 = mz × sqrt(2gH1) / G  =>  H1 = (v4 × G / mz)² / (2g)
H1 = (v4 * G / mz)^2 / (2 * g);

fprintf('\n');
fprintf('式(17): H1 = (v4×G/mz)² / (2g)\n');
fprintf('      = (%.4f × %.4f / %.0f)² / (2 × %.1f)\n', v4, G, mz, g);
fprintf('>>> 斜抛运动高度 H1 = %.4f mm\n', H1 * 1000);
fprintf('\n');

%% ========== 步骤10: 计算跳跃总高度 H (式4) ==========
fprintf('>>> 步骤10: 计算跳跃总高度 H (式4)\n');
fprintf('--------------------------------------------\n');

% 式(4): H = H1 + H2
H = H1 + H2;

fprintf('H = H1 + H2\n');
fprintf('  = %.4f mm + %.4f mm\n', H1 * 1000, H2 * 1000);
fprintf('>>> 跳跃总高度 H = %.4f mm\n', H * 1000);
fprintf('\n');

%% ========== 步骤11: 验证角速度约束 ==========
fprintf('>>> 步骤11: 验证角速度约束 (式32, 33a, 33b)\n');
fprintf('--------------------------------------------\n');

% 式(32): ω0 = (ε1-ε2)×F̄ / (mz×√(2gH1))
% 式(33a): ω1 = (ε3-ε2)×√(g/(2H1))
% 式(33b): ω2 = (δ3-δ2)×√(g/(2H1))

omega0 = abs((eps1 - eps2) * F_bar / (mz * sqrt(2 * g * H1)));
omega1 = abs((eps3 - eps2) * sqrt(g / (2 * H1)));
omega2 = abs((delta3 - delta2) * sqrt(g / (2 * H1)));

fprintf('ω0 = |ε1-ε2|×F̄ / (mz×√(2gH1)) = %.4f rad/s\n', omega0);
fprintf('ω1 = |ε3-ε2|×√(g/(2H1)) = %.4f rad/s\n', omega1);
fprintf('ω2 = |δ3-δ2|×√(g/(2H1)) = %.4f rad/s\n', omega2);
fprintf('最大允许转速 ωmax = %.0f rad/s\n', w_max);
fprintf('\n');

if omega0 <= w_max && omega1 <= w_max && omega2 <= w_max
    fprintf('>>> 所有角速度约束满足 ✓\n');
else
    fprintf('>>> 角速度约束不满足 ✗\n');
end
fprintf('\n');

%% ========== 步骤12: 约束条件检查 ==========
fprintf('>>> 步骤12: 约束条件检查 (仅显示，不影响计算)\n');
fprintf('--------------------------------------------\n');

% 安全系数 (与优化代码一致)
ka13 = 1.3;  % 地面缩腿和飞行阶段的膝关节角度安全系数
ka2 = 2;     % 刚离地阶段的膝关节角度安全系数

% --- 约束1: k范围约束 ---
k_lower = 0.2; k_upper = 0.8;
fprintf('约束1 - k范围: %.2f < k < %.2f\n', k_lower, k_upper);
fprintf('  当前值: k = %.4f\n', k);
if k > k_lower && k < k_upper
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束2: ε1范围约束 ---
% 计算角度限制 ε_13lim (式35a)
sum_r = r1 + r2;
if sum_r + l1 <= l2
    eps_13lim = asin(r1 / l1);
elseif sum_r + l2 <= l1
    eps_13lim = asin(r2 / l2);
else
    eps_13lim = acos((l1^2 + l2^2 - sum_r^2) / (2*l1*l2));
end
eps1_lower = pi/2;
eps1_upper = pi - eps_13lim * ka13;
fprintf('\n约束2 - ε1范围: %.2f° < ε1 < %.2f°\n', rad2deg(eps1_lower), rad2deg(eps1_upper));
fprintf('  当前值: ε1 = %.4f°\n', eps1_deg);
if eps1 > eps1_lower && eps1 < eps1_upper
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束3: ε2范围约束 ---
eps_2lim = asin(r1 / l1);
eps2_lower = eps_2lim * ka2;
eps2_upper = pi/2;
fprintf('\n约束3 - ε2范围: %.2f° < ε2 < %.2f°\n', rad2deg(eps2_lower), rad2deg(eps2_upper));
fprintf('  当前值: ε2 = %.4f°\n', eps2_deg);
if eps2 > eps2_lower && eps2 < eps2_upper
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束4: ε3范围约束 ---
fprintf('\n约束4 - ε3范围: ε1 ≤ ε3 < %.2f°\n', rad2deg(eps1_upper));
fprintf('  当前值: ε3 = %.4f°, ε1 = %.4f°\n', eps3_deg, eps1_deg);
if eps3 >= eps1 && eps3 < eps1_upper
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束5: δ3范围约束 ---
% 计算 δ_3lim (式35b)
if sum_r + l1 <= l2
    delta_3lim = 0;
else
    delta_3lim = pi/2 - eps_13lim * ka13;
end
delta3_upper = pi/2 + delta_3lim;
fprintf('\n约束5 - δ3范围: δ2 < δ3 ≤ %.2f°\n', rad2deg(delta3_upper));
fprintf('  当前值: δ3 = %.4f°, δ2 = %.4f°\n', delta3_deg, delta2_deg);
if delta3 > delta2 && delta3 <= delta3_upper
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束6: ε1 ≤ π/2 + δ1 ---
eps1_limit = pi/2 + delta1;
fprintf('\n约束6 - ε1姿态约束: ε1 ≤ π/2 + δ1 = %.2f°\n', rad2deg(eps1_limit));
fprintf('  当前值: ε1 = %.4f°\n', eps1_deg);
if eps1 <= eps1_limit
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束7: ε3 ≤ π/2 + δ3 ---
eps3_limit = pi/2 + delta3;
fprintf('\n约束7 - ε3姿态约束: ε3 ≤ π/2 + δ3 = %.2f°\n', rad2deg(eps3_limit));
fprintf('  当前值: ε3 = %.4f°\n', eps3_deg);
if eps3 <= eps3_limit
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束8: H2范围约束 ---
H2_upper = l1 + l2 - r2;
fprintf('\n约束8 - H2范围: 0 < H2 ≤ %.2f mm\n', H2_upper*1000);
fprintf('  当前值: H2 = %.4f mm\n', H2*1000);
if H2 > 0 && H2 <= H2_upper
    fprintf('  状态: ✓ 满足\n');
else
    fprintf('  状态: ✗ 超出限度!\n');
end

% --- 约束9: F̄范围约束 ---
fprintf('\n约束9 - F̄范围: 0 < F̄ ≤ %.2f N\n', F_upper);
fprintf('  理论值: F̄_theory = %.2f N\n', F_bar_theory);
fprintf('  实际使用值: F̄ = %.2f N (取min后)\n', F_bar);
if F_bar_theory > 0 && F_bar_theory <= F_upper
    fprintf('  理论值状态: ✓ 满足\n');
else
    fprintf('  理论值状态: ✗ 超出限度! (理论值%.2f > 上限%.2f)\n', F_bar_theory, F_upper);
end

% --- 约束10: 角速度约束 ---
fprintf('\n约束10 - 角速度约束: |ω| ≤ %.0f rad/s\n', w_max);
fprintf('  ω0 = %.4f rad/s', omega0);
if omega0 <= w_max
    fprintf(' ✓\n');
else
    fprintf(' ✗ 超出!\n');
end
fprintf('  ω1 = %.4f rad/s', omega1);
if omega1 <= w_max
    fprintf(' ✓\n');
else
    fprintf(' ✗ 超出!\n');
end
fprintf('  ω2 = %.4f rad/s', omega2);
if omega2 <= w_max
    fprintf(' ✓\n');
else
    fprintf(' ✗ 超出!\n');
end

fprintf('\n');

%% ========== 结果汇总与验证 ==========
fprintf('============================================\n');
fprintf('             计算结果汇总\n');
fprintf('============================================\n');
fprintf('\n');
fprintf('+---------------------------+-------------+-------------+\n');
fprintf('|         参数名称          |  计算值     |  文档值     |\n');
fprintf('+---------------------------+-------------+-------------+\n');
fprintf('| 空中主动缩腿高度 H2       | %8.1f mm | %8.1f mm |\n', H2*1000, 259.9);
fprintf('| 斜抛运动高度 H1           | %8.1f mm | %8.1f mm |\n', H1*1000, 280.5-259.9);
fprintf('| 跳跃总高度 H = H1+H2      | %8.1f mm | %8.1f mm |\n', H*1000, 280.5);
fprintf('| 地面伸腿高度(重心上升)    | %8.1f mm | %8.1f mm |\n', h_ground_extend*1000, 222.1);
fprintf('+---------------------------+-------------+-------------+\n');
fprintf('\n');

% 计算误差
error_H2 = abs(H2*1000 - 259.9);
error_H1 = abs(H1*1000 - (280.5-259.9));
error_H = abs(H*1000 - 280.5);
error_h_ground = abs(h_ground_extend*1000 - 222.1);

fprintf('与文档值的误差:\n');
fprintf('  H2误差: %.2f mm (%.2f%%)\n', error_H2, error_H2/259.9*100);
fprintf('  H1误差: %.2f mm (%.2f%%)\n', error_H1, error_H1/20.6*100);
fprintf('  H 误差: %.2f mm (%.2f%%)\n', error_H, error_H/280.5*100);
fprintf('  地面伸腿高度误差: %.2f mm (%.2f%%)\n', error_h_ground, error_h_ground/222.1*100);
fprintf('\n');

fprintf('说明:\n');
fprintf('  - H2(空中主动缩腿高度): 空中收腿增加的跳跃高度\n');
fprintf('  - H1(斜抛运动高度): 质心抛射运动产生的高度\n');
fprintf('  - H(跳跃总高度) = H1 + H2\n');
fprintf('  - 地面伸腿高度: 起跳过程中重心上升的几何高度(非H1)\n');
fprintf('\n');

% 最终验证结论
fprintf('============================================\n');
fprintf('             验证结论\n');
fprintf('============================================\n');
if error_H < 5  % 允许5mm误差
    fprintf('>>> 验证成功! 计算结果与文档给出的280.5mm基本一致\n');
else
    fprintf('>>> 验证存在偏差，可能需要检查公式理解或参数设置\n');
end
fprintf('============================================\n');

