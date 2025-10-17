%% 平均输出扭矩 t 表达式


%% 计算简化的雅可比
k_rky = diff(rky_sol, theta2); % 雅可比矩阵的 j22
tau = m * g / 2 * k_rky / 1000 / 2; % 求关节扭矩

b_ang = -30;
theta2_val = deg2rad(b_ang);

tau_val = subs(tau, [l1, l2], [l1_val, l2_val]);
t_ave = [];
while b_ang <= 60
    theta2_val = deg2rad(b_ang);
    tau_num = double(subs(tau_val, theta2, theta2_val));
    t_ave(end + 1) = abs(tau_num);
    b_ang = b_ang + 2;
end
t_average(i) = mean(t_ave);

%% 计算行程中 t 的极值和所在角度
j_num = subs(J, [l1, l2, k], [l1_val, l2_val, k_val]);
interval = 5;
a_ang = 60;
b_ang = -30;
t = [];
while a_ang > 15
    while b_ang < 15
        theta1_val = deg2rad(180 - a_ang);
        theta2_val = deg2rad(b_ang);
        [theta1, theta2] = [theta1_val, theta2_val];
        J_final = double(subs(j_num, [theta1, theta2], [theta1_val, theta2_val]));
        F = [0; m * g / 2]; % 地面反作用力(N)
        tau_final = J_final' * F; % 关节扭矩(N·m)
        t(end + 1) = abs(tau_final(1));
        t(end + 1) = abs(tau_final(2));
        b_ang = b_ang + interval;
    end
    b_ang = -20;
    a_ang = a_ang - interval;
end

[max_t, index_max] = max(t);
result_max = findindex(index_max);
fprintf('最大值: %.3f, A 关节角: %.f, B 关节角: %.f。\n', max_t, result_max(1), result_max(2));

[min_t, index_min] = min(t);
result_min = findindex(index_min);
fprintf('最小值: %.3f, A 关节角: %.f, B 关节角: %.f。\n', min_t, result_min(1), result_min(2));

%% findindex 函数
function result = findindex(index)
    if rem(index, 2) == 0
        index = index / 2;
    else
        index = (index + 1) / 2;
    end
    index_b = mod(index, 10); % 取余
    index_a = (index - index_b) / 10; % 取商
    aang = 80 - 5 * index_a;
    bang = -20 + 5 * (index_b - 1);
    result = [aang, bang];
end