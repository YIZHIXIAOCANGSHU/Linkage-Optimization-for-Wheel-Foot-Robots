%% 两轮足机器人跳跃高度PSO优化
% 基于专利 CN 116738739 A 的数学模型
% 完整的从计算到优化的参数寻优代码

clear; clc; close all;

%% ========== 步骤1: 机器人参数初始化 (表1) ==========
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

% 安全系数
ka13 = 1.3;         % 地面缩腿和飞行阶段的膝关节角度安全系数
ka2 = 2;            % 刚离地阶段的膝关节角度安全系数

% 物理常数
g = 9.8;            % 重力加速度 (m/s^2)

% 将参数打包为结构体，便于传递
params.m1 = m1; params.m2 = m2; params.m3 = m3; params.m4 = m4; params.mz = mz;
params.lz = lz; params.r1 = r1; params.r2 = r2;
params.b = b; params.c = c;
params.M2 = M2; params.w_max = w_max;
params.ka13 = ka13; params.ka2 = ka2;
params.g = g;

%% ========== 步骤2: 定义辅助函数 ==========

% 式(2): 计算髋关节角度δ与膝关节角度ε的关系
% δ = arccot[csc(ε)(cos(ε) + (m3*c+m4)*l2 / ((m2*b+m3+m4)*l1))]
calc_delta = @(eps, k) calc_delta_func(eps, k, params);

% 计算角度限制 ε_13lim (式35a)
calc_eps_13lim = @(k) calc_eps_13lim_func(k, params);

% 计算角度限制 δ_3lim (式35b)
calc_delta_3lim = @(k, eps_13lim) calc_delta_3lim_func(k, eps_13lim, params);

%% ========== 步骤3: PSO优化设置 ==========

% 优化变量: x = [k, eps1, eps2, eps3, delta3]
% k: 小腿长度占比
% eps1: 地面缩腿阶段膝关节角度 (rad)
% eps2: 刚离地阶段膝关节角度 (rad)
% eps3: 飞行阶段膝关节角度 (rad)
% delta3: 飞行阶段髋关节角度 (rad)

% 变量数量
nvars = 5;

% 变量边界 (宽松边界，实际约束在目标函数中处理)
lb = [0.2, pi/2, 0.1, pi/2, 0.1];
ub = [0.8, pi, pi/2, pi, pi];

% 目标函数 (最大化跳跃高度 -> 最小化负的跳跃高度)
objective = @(x) objective_function(x, params);

%% ========== 步骤4: 多策略PSO优化（减少局部最优） ==========
fprintf('========== 开始多策略PSO优化 ==========\n');
fprintf('优化变量: [k, ε1, ε2, ε3, δ3]\n');
fprintf('目标: 最大化跳跃总高度 H\n');
fprintf('F̄约束上限: (mz-m1)×g = %.2f N\n', (mz-m1)*g);
fprintf('策略: 多次独立运行 + 混沌初始化 + 自适应参数\n\n');

% ===== 优化参数设置 =====
num_runs = 5;           % 独立运行次数
swarm_size = 300;       % 每次运行的粒子数量
max_iter = 500;         % 每次运行的最大迭代次数

% 用于记录所有运行结果
all_results = zeros(num_runs, nvars);
all_fvals = zeros(num_runs, 1);
all_histories = cell(num_runs, 1);

% 全局变量用于记录收敛历史
global current_run_history;

% 记录总开始时间
tic;

% 自适应PSO参数（根据运行次数调整）
% 不同运行使用不同的参数组合，增加多样性
inertia_ranges = {[0.4, 0.9], [0.3, 0.8], [0.5, 1.0], [0.2, 0.7], [0.4, 0.95]};
c1_values = [1.49, 1.7, 1.3, 1.5, 1.6];
c2_values = [1.49, 1.3, 1.7, 1.6, 1.4];

for run_idx = 1:num_runs
    fprintf('>>> 第 %d/%d 次运行 (惯性:[%.1f,%.1f], c1=%.2f, c2=%.2f)...\n', ...
        run_idx, num_runs, inertia_ranges{run_idx}(1), inertia_ranges{run_idx}(2), ...
        c1_values(run_idx), c2_values(run_idx));
    
    % 生成混沌初始种群（使用Logistic映射）
    initial_swarm = generate_chaotic_swarm(swarm_size, nvars, lb, ub, run_idx);
    
    % 添加文档最优参数作为初始粒子之一（第一次运行）
    if run_idx == 1
        doc_params = [0.4675, deg2rad(147.3934), deg2rad(69.4081), ...
                      deg2rad(147.3934), deg2rad(114.8207)];
        initial_swarm(1, :) = doc_params;
    end
    
    % 重置收敛历史
    current_run_history = [];
    
    options = optimoptions('particleswarm', ...
        'SwarmSize', swarm_size, ...
        'MaxIterations', max_iter, ...
        'MaxStallIterations', 80, ...
        'FunctionTolerance', 1e-12, ...
        'InertiaRange', inertia_ranges{run_idx}, ...
        'SelfAdjustmentWeight', c1_values(run_idx), ...
        'SocialAdjustmentWeight', c2_values(run_idx), ...
        'MinNeighborsFraction', 0.15 + 0.05*run_idx, ...
        'InitialSwarmMatrix', initial_swarm, ...
        'Display', 'off', ...
        'OutputFcn', @pso_history_func);
    
    % 执行PSO
    [x_run, fval_run, ~, ~] = particleswarm(objective, nvars, lb, ub, options);
    
    % 记录结果
    all_results(run_idx, :) = x_run;
    all_fvals(run_idx) = fval_run;
    all_histories{run_idx} = current_run_history;
    
    fprintf('   完成! H = %.2f mm\n', -fval_run * 1000);
end

elapsed_time = toc;

% 找到最优结果
[best_fval, best_idx] = min(all_fvals);
x_opt = all_results(best_idx, :);
fval = best_fval;

% 合并收敛历史（使用最优运行的历史）
convergence_history = all_histories{best_idx};

fprintf('\n========== 多次运行结果汇总 ==========\n');
for i = 1:num_runs
    marker = '';
    if i == best_idx
        marker = ' <-- 最优';
    end
    fprintf('运行 %d: H = %.2f mm%s\n', i, -all_fvals(i) * 1000, marker);
end
fprintf('结果标准差: %.2f mm\n', std(-all_fvals * 1000));

%% ========== 步骤5: 输出优化结果 ==========
fprintf('\n========== 最优结果 (来自运行 %d) ==========\n', best_idx);
fprintf('总优化耗时: %.2f 秒\n', elapsed_time);
fprintf('独立运行次数: %d\n', num_runs);
fprintf('每次粒子数: %d, 最大迭代: %d\n', swarm_size, max_iter);

% 提取最优参数
k_opt = x_opt(1);
eps1_opt = x_opt(2);
eps2_opt = x_opt(3);
eps3_opt = x_opt(4);
delta3_opt = x_opt(5);

% 计算对应的髋关节角度
l1_opt = k_opt * lz;
l2_opt = (1 - k_opt) * lz;
delta1_opt = calc_delta(eps1_opt, k_opt);
delta2_opt = calc_delta(eps2_opt, k_opt);

% 计算跳跃高度各分量
[H_total, H1, H2, F_bar] = calc_jump_height_detail(x_opt, params);

fprintf('\n--- 最优参数 ---\n');
fprintf('小腿长度占比 k = %.4f\n', k_opt);
fprintf('小腿长度 l1 = %.2f mm\n', l1_opt * 1000);
fprintf('大腿长度 l2 = %.2f mm\n', l2_opt * 1000);
fprintf('\n--- 角度参数 (度) ---\n');
fprintf('地面缩腿膝关节角度 ε1 = %.4f°\n', rad2deg(eps1_opt));
fprintf('刚离地膝关节角度 ε2 = %.4f°\n', rad2deg(eps2_opt));
fprintf('飞行阶段膝关节角度 ε3 = %.4f°\n', rad2deg(eps3_opt));
fprintf('地面缩腿髋关节角度 δ1 = %.4f°\n', rad2deg(delta1_opt));
fprintf('刚离地髋关节角度 δ2 = %.4f°\n', rad2deg(delta2_opt));
fprintf('飞行阶段髋关节角度 δ3 = %.4f°\n', rad2deg(delta3_opt));
fprintf('\n--- 跳跃高度 (mm) ---\n');
fprintf('斜抛运动高度 H1 = %.1f mm\n', H1 * 1000);
fprintf('主动缩腿高度 H2 = %.1f mm\n', H2 * 1000);
fprintf('跳跃总高度 H = %.1f mm\n', H_total * 1000);
fprintf('\n--- 其他参数 ---\n');
fprintf('平均作用力 F̄ = %.2f N (取min后的实际值)\n', F_bar);
fprintf('F̄约束上限 = %.2f N\n', (mz-m1)*g);

%% ========== 步骤6: 绘制PSO收敛曲线 ==========
fprintf('\n========== 绘制收敛曲线 ==========\n');

figure('Name', 'PSO多策略优化收敛曲线', 'NumberTitle', 'off', 'Position', [100, 100, 1000, 700]);

% 子图1：所有运行的收敛曲线对比
subplot(2, 2, [1, 2]);
colors = lines(num_runs);
legend_entries = cell(num_runs, 1);
hold on;
for i = 1:num_runs
    if ~isempty(all_histories{i})
        h = plot(all_histories{i}(:, 1), all_histories{i}(:, 2), '-', 'LineWidth', 1.5, 'Color', colors(i, :));
        if i == best_idx
            set(h, 'LineWidth', 2.5);
            legend_entries{i} = sprintf('运行%d: %.1fmm (最优)', i, -all_fvals(i)*1000);
        else
            legend_entries{i} = sprintf('运行%d: %.1fmm', i, -all_fvals(i)*1000);
        end
    end
end
hold off;
xlabel('迭代次数', 'FontSize', 11);
ylabel('最优跳跃高度 H (mm)', 'FontSize', 11);
title('多次独立运行收敛曲线对比', 'FontSize', 13);
legend(legend_entries, 'Location', 'southeast', 'FontSize', 9);
grid on;

% 子图2：最优运行的收敛曲线
subplot(2, 2, 3);
if ~isempty(convergence_history)
    plot(convergence_history(:, 1), convergence_history(:, 2), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(convergence_history(end, 1), convergence_history(end, 2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    hold off;
    xlabel('迭代次数', 'FontSize', 11);
    ylabel('最优跳跃高度 H (mm)', 'FontSize', 11);
    title(sprintf('最优运行 (运行%d) 收敛曲线', best_idx), 'FontSize', 12);
    legend({'收敛曲线', sprintf('最优: %.2fmm', H_total*1000)}, 'Location', 'southeast');
    grid on;
end

% 子图3：各运行最终结果柱状图
subplot(2, 2, 4);
bar_colors = repmat([0.3, 0.6, 0.9], num_runs, 1);
bar_colors(best_idx, :) = [0.9, 0.3, 0.3];  % 最优用红色
b = bar(-all_fvals * 1000);
b.FaceColor = 'flat';
b.CData = bar_colors;
xlabel('运行次数', 'FontSize', 11);
ylabel('跳跃高度 H (mm)', 'FontSize', 11);
title('各运行最终结果对比', 'FontSize', 12);
% 添加数值标签
for i = 1:num_runs
    text(i, -all_fvals(i)*1000 + 5, sprintf('%.1f', -all_fvals(i)*1000), ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end
grid on;

% 添加总标题
sgtitle(sprintf('两轮足机器人跳跃高度 - 多策略PSO优化\n最优高度: %.2f mm | 耗时: %.2f s | 运行次数: %d', ...
    H_total * 1000, elapsed_time, num_runs), 'FontSize', 14);

fprintf('收敛曲线已绘制完成\n');
if ~isempty(convergence_history)
    fprintf('最优运行初始值: %.2f mm\n', convergence_history(1, 2));
    fprintf('最优运行最终值: %.2f mm\n', convergence_history(end, 2));
    if convergence_history(1, 2) > 0
        fprintf('最优运行提升: %.2f mm (%.2f%%)\n', ...
            convergence_history(end, 2) - convergence_history(1, 2), ...
            (convergence_history(end, 2) - convergence_history(1, 2)) / convergence_history(1, 2) * 100);
    end
end

%% ========== 辅助函数定义 ==========

function stop = pso_history_func(optimValues, state)
    % PSO输出函数：记录每次迭代的最优值
    global current_run_history;
    stop = false;
    if strcmp(state, 'iter')
        current_run_history = [current_run_history; ...
            optimValues.iteration, -optimValues.bestfval * 1000];
    end
end

function swarm = generate_chaotic_swarm(pop_size, dim, lb, ub, seed)
    % 使用Logistic映射生成混沌初始种群
    rng(seed * 1000);  % 设置随机种子保证可重复性
    
    % 初始化混沌序列
    chaos = zeros(pop_size, dim);
    x0 = 0.1 + 0.8 * rand(1, dim);  % 随机初始值
    
    % Logistic映射参数
    mu = 4;  % 混沌参数
    
    for d = 1:dim
        x = x0(d);
        for i = 1:pop_size
            x = mu * x * (1 - x);  % Logistic映射
            chaos(i, d) = x;
        end
    end
    
    % 映射到搜索空间
    swarm = lb + chaos .* (ub - lb);
    
    % 添加一些均匀分布的粒子增加覆盖性
    num_uniform = round(pop_size * 0.2);
    uniform_particles = lb + rand(num_uniform, dim) .* (ub - lb);
    swarm(end-num_uniform+1:end, :) = uniform_particles;
end

function delta = calc_delta_func(eps, k, params)
    % 式(2): 计算髋关节角度δ
    % δ = arccot[csc(ε)(cos(ε) + (m3*c+m4)*l2 / ((m2*b+m3+m4)*l1))]
    
    l1 = k * params.lz;
    l2 = (1 - k) * params.lz;
    
    m2 = params.m2; m3 = params.m3; m4 = params.m4;
    b = params.b; c = params.c;
    
    % 计算括号内的表达式
    term = csc(eps) * (cos(eps) + (m3*c + m4)*l2 / ((m2*b + m3 + m4)*l1));
    
    % arccot(x) = atan(1/x)
    delta = atan(1 / term);
    
    % 确保角度在正确范围内
    if delta < 0
        delta = delta + pi;
    end
end

function eps_13lim = calc_eps_13lim_func(k, params)
    % 式(35a): 计算膝关节角度限制 ε_13lim
    l1 = k * params.lz;
    l2 = (1 - k) * params.lz;
    r1 = params.r1;
    r2 = params.r2;
    
    sum_r = r1 + r2;
    
    if sum_r + l1 <= l2
        eps_13lim = asin(r1 / l1);
    elseif sum_r + l2 <= l1
        eps_13lim = asin(r2 / l2);
    else
        eps_13lim = acos((l1^2 + l2^2 - sum_r^2) / (2*l1*l2));
    end
end

function delta_3lim = calc_delta_3lim_func(k, eps_13lim, params)
    % 式(35b): 计算髋关节角度限制 δ_3lim
    l1 = k * params.lz;
    l2 = (1 - k) * params.lz;
    r1 = params.r1;
    r2 = params.r2;
    ka13 = params.ka13;
    
    sum_r = r1 + r2;
    
    if sum_r + l1 <= l2
        delta_3lim = 0;
    elseif sum_r + l2 <= l1
        delta_3lim = pi/2 - eps_13lim * ka13;
    else
        delta_3lim = pi/2 - eps_13lim * ka13;
    end
end

function f = objective_function(x, params)
    % 目标函数: 返回负的跳跃高度 (因为PSO是最小化)
    % 如果违反约束，返回一个很大的惩罚值
    
    % 提取变量
    k = x(1);
    eps1 = x(2);
    eps2 = x(3);
    eps3 = x(4);
    delta3 = x(5);
    
    % 计算腿长
    l1 = k * params.lz;
    l2 = (1 - k) * params.lz;
    
    % 计算髋关节角度
    delta1 = calc_delta_func(eps1, k, params);
    delta2 = calc_delta_func(eps2, k, params);
    
    % ===== 检查约束条件 (式34a, 34b, 35a, 35b) =====
    
    % 计算角度限制
    eps_13lim = calc_eps_13lim_func(k, params);
    delta_3lim = calc_delta_3lim_func(k, eps_13lim, params);
    
    % 约束1: 0.2 < k < 0.8
    if k <= 0.2 || k >= 0.8
        f = 1e6;
        return;
    end
    
    % 约束2: π/2 < ε1 < π - ε_13lim * ka13
    eps1_lower = pi/2;
    eps1_upper = pi - eps_13lim * params.ka13;
    if eps1 <= eps1_lower || eps1 >= eps1_upper
        f = 1e6;
        return;
    end
    
    % 约束3: ε_2lim * ka2 < ε2 < π/2
    % ε_2lim = arcsin(r1/l1) 刚离地时膝关节的机械限制角度
    eps_2lim = asin(params.r1 / l1);
    eps2_lower = eps_2lim * params.ka2;
    eps2_upper = pi/2;
    if eps2 <= eps2_lower || eps2 >= eps2_upper
        f = 1e6;
        return;
    end
    
    % 约束4: ε1 ≤ ε3 < π - ε_13lim * ka13
    if eps3 < eps1 || eps3 >= eps1_upper
        f = 1e6;
        return;
    end
    
    % 约束5: δ2 < δ3 ≤ π/2 + δ_3lim
    delta3_upper = pi/2 + delta_3lim;
    if delta3 <= delta2 || delta3 > delta3_upper
        f = 1e6;
        return;
    end
    
    % 约束6: ε1 ≤ π/2 + δ1
    if eps1 > pi/2 + delta1
        f = 1e6;
        return;
    end
    
    % 约束7: ε3 ≤ π/2 + δ3
    if eps3 > pi/2 + delta3
        f = 1e6;
        return;
    end
    
    % ===== 计算跳跃高度 =====
    [H, H1, H2, F_bar] = calc_jump_height_detail(x, params);
    
    % 检查计算结果是否有效
    if ~isreal(H) || isnan(H) || isinf(H) || H <= 0
        f = 1e6;
        return;
    end
    
    % 约束8: 0 < H2 ≤ l1 + l2 - r2
    H2_upper = l1 + l2 - params.r2;
    if H2 <= 0 || H2 > H2_upper
        f = 1e6;
        return;
    end
    
    % 约束9: F̄ > 0 (F̄已在calc_jump_height_detail中取min限制上限)
    if F_bar <= 0
        f = 1e6;
        return;
    end
    
    % 约束10: 角速度约束 |ω0|, |ω1|, |ω2| ≤ ω_max
    % 式(32): ω0 = (ε1-ε2)*F̄ / (mz*√(2gH1))
    % 式(33a): ω1 = (ε3-ε2)*√(g/(2H1))
    % 式(33b): ω2 = (δ3-δ2)*√(g/(2H1))
    if H1 > 0
        omega0 = abs((eps1 - eps2) * F_bar / (params.mz * sqrt(2 * params.g * H1)));
        omega1 = abs((eps3 - eps2) * sqrt(params.g / (2 * H1)));
        omega2 = abs((delta3 - delta2) * sqrt(params.g / (2 * H1)));
        
        if omega0 > params.w_max || omega1 > params.w_max || omega2 > params.w_max
            f = 1e6;
            return;
        end
    end
    
    % 返回负的跳跃高度 (最小化负值 = 最大化正值)
    f = -H;
end

function [H, H1, H2, F_bar] = calc_jump_height_detail(x, params)
    % 计算跳跃高度的详细函数 (与jump_single_calculation.m保持一致)
    % 返回: H(总高度), H1(斜抛高度), H2(主动缩腿高度), F_bar(平均作用力)
    
    % 提取变量
    k = x(1);
    eps1 = x(2);
    eps2 = x(3);
    eps3 = x(4);
    delta3 = x(5);
    
    % 提取参数
    m1 = params.m1; m2 = params.m2; m3 = params.m3; m4 = params.m4; mz = params.mz;
    lz = params.lz; b = params.b; c = params.c;
    M2 = params.M2; g = params.g;
    
    % 计算腿长
    l1 = k * lz;
    l2 = (1 - k) * lz;
    
    % 计算髋关节角度 (式2)
    delta1 = calc_delta_func(eps1, k, params);
    delta2 = calc_delta_func(eps2, k, params);
    
    % ===== 步骤1: 计算H2 - 主动缩腿高度 (式5) =====
    % H2 = l1[cos(ε2-δ2) - cos(ε3-δ3)] + l2(cosδ2 - cosδ3)
    term_H2_1 = l1 * (cos(eps2 - delta2) - cos(eps3 - delta3));
    term_H2_2 = l2 * (cos(delta2) - cos(delta3));
    H2 = term_H2_1 + term_H2_2;
    
    % ===== 步骤2: 计算平均作用力 F̄ =====
    % 式(22)计算理论值，约束(23)限制上限: 0 < F̄ ≤ (mz-m1)*g
    sin_diff1 = sin(eps1 - delta1);
    sin_diff2 = sin(eps2 - delta2);
    
    if abs(sin_diff2 - sin_diff1) < 1e-10 || sin_diff1 <= 0 || sin_diff2 <= 0
        H = 0; H1 = 0; F_bar = 0;
        return;
    end
    
    % 式(22)计算理论F̄值
    F_bar_theory = 2 * M2 / (l1 * (sin_diff2 - sin_diff1)) * log(sin_diff2 / sin_diff1);
    
    % 约束上限 (式23)
    F_upper = (mz - m1) * g;
    
    % 实际F̄取理论值和约束上限的较小值
    F_bar = min(F_bar_theory, F_upper);
    
    % 如果F_bar为负或为0，返回无效值
    if F_bar <= 0
        H = 0; H1 = 0;
        return;
    end
    
    % ===== 步骤3: 计算G系数 (式18) =====
    % G = (m2*b+m3)*cos²δ2 + 0.5*(m3*(c-1)-m2*b)*cot(ε2)*sin(2δ2) + m3*c*sin²δ2 + m4
    cos_d2 = cos(delta2);
    sin_d2 = sin(delta2);
    cot_e2 = cot(eps2);
    sin_2d2 = sin(2 * delta2);
    
    term_G1 = (m2*b + m3) * cos_d2^2;
    term_G2 = 0.5 * (m3*(c-1) - m2*b) * cot_e2 * sin_2d2;
    term_G3 = m3 * c * sin_d2^2;
    term_G4 = m4;
    G = term_G1 + term_G2 + term_G3 + term_G4;
    
    % ===== 步骤4: 计算K系数 (式25) =====
    % K = m2*b²*(cos²δ2/sin²ε2) + m3*cos²δ2*[1+(c*tanδ2-(1-c)*cotε2)²] + m4
    sin_e2 = sin(eps2);
    tan_d2 = tan(delta2);
    
    term_K1 = m2 * b^2 * (cos_d2^2 / sin_e2^2);
    inner_term = c * tan_d2 - (1 - c) * cot_e2;
    term_K2 = m3 * cos_d2^2 * (1 + inner_term^2);
    term_K3 = m4;
    K = term_K1 + term_K2 + term_K3;
    
    % ===== 步骤5: 计算T和Q (式25) =====
    % T = cosδ2 - cosδ1
    % Q = cos(ε2-δ2) - cos(ε1-δ1)
    T = cos(delta2) - cos(delta1);
    Q = cos(eps2 - delta2) - cos(eps1 - delta1);
    
    % ===== 步骤6: 计算斜抛运动高度H1 (式24 + 式17) =====
    % 式(24): v4² = (2/K) × {[F̄-(m2*b+m3+m4)*g]*l1*Q + [F̄-(m3*c+m4)*g]*l2*T}
    bracket1 = F_bar - (m2*b + m3 + m4) * g;
    bracket2 = F_bar - (m3*c + m4) * g;
    
    term1 = bracket1 * l1 * Q;
    term2 = bracket2 * l2 * T;
    v4_squared = (2/K) * (term1 + term2);
    
    if v4_squared <= 0
        H = 0; H1 = 0;
        return;
    end
    
    v4 = sqrt(v4_squared);
    
    % 式(17): v4 = mz × sqrt(2gH1) / G  =>  H1 = (v4 × G / mz)² / (2g)
    H1 = (v4 * G / mz)^2 / (2 * g);
    
    % ===== 步骤7: 计算跳跃总高度 (式4) =====
    % H = H1 + H2
    H = H1 + H2;
end

