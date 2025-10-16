%% 四连杆末端轨迹优化 - fmincon SQP 算法
% 四连杆机构尺寸优化，只使用 fmincon SQP 算法
clc;
clear;

% 定义机构的初始杆长
l4 = 0.077*sqrt(2); %初始值
l2 = 0.22;
l3 = 0.223;
l23 = 0.06;
l1 = 0.23;

initParams = [l1,l2,l3,l4,l23];

% 定义优化约束，确保杆长在合理范围内
lb = [0.16, 0.2, 0.2, 0.105, 0.05]; % 下限
ub = [0.28, 0.28, 0.28, 0.12, 0.08]; % 上限

% 初始化优化过程记录变量
global optHistory iterationCount lastPlotTime
optHistory = [];
iterationCount = 0;
lastPlotTime = tic;

% 创建图形窗口
figure(1);
set(1, 'Name', '优化过程可视化', 'NumberTitle', 'off');

% 使用 fmincon SQP 算法
fprintf('正在运行 fmincon SQP 算法...\n');
options_fmincon = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', ...
    'PlotFcn', @animateOptimization);
[optimizedParams_fmincon, fval_fmincon] = fmincon(@objectiveFun, initParams, [], [], [], [], lb, ub, [], options_fmincon);

fprintf('fmincon 结果: 目标函数值 = %.6f\n', fval_fmincon);
fprintf('优化参数: %s\n', mat2str(optimizedParams_fmincon, 4));

%% 优化目标函数：末端轨迹尽可能接近直线
function error = objectiveFun(params)
    % 检查参数是否有效
    if any(isnan(params)) || any(isinf(params)) || any(params <= 0)
        error = inf;
        return;
    end
    
    L1 = params(1);
    L2 = params(2);
    L3 = params(3);
    L4 = params(4);
    L23 = params(5);

    % 定义输入角度范围
    theta = linspace(pi/6, pi/3, 100);

    % 计算末端轨迹
    endEffectorX = zeros(size(theta));

    for i = 1:length(theta)
        % 计算末端位置误差，此处的计算需根据具体四连杆机构的运动学方程设置
        endEffectorX(i) = End_Pos_calc(L1, L2, L3, L4, L23, theta(i));
    end

    % 计算误差平方和作为目标函数，确保返回实数值
    error = sum(real(endEffectorX).^2);
    
    % 确保返回的值是实数且非负
    if isnan(error) || isinf(error)
        error = realmax;
    end
    error = max(0, real(error));
end

% 末端位置计算函数（需根据实际四连杆结构设计）
function [Ex] = End_Pos_calc(L1, L2, L3, L4, L23, theta)
    try
        % 实际结构需要根据四连杆几何关系进行计算
        % 四连杆结构的简单示例计算
        angBAD = pi/4 + theta;
        
        % 确保cos参数在有效范围内
        cos_angBAD = cos(angBAD);
        cos_angBAD(cos_angBAD > 1) = 1;
        cos_angBAD(cos_angBAD < -1) = -1;
        
        BD = sqrt(L2.^2+L4.^2-2.*L2.*L4.*cos_angBAD);
        
        % 检查BD是否有效
        BD(isnan(BD) | isinf(BD)) = 0;
        BD(BD < 0) = 0;
        
        % 计算余弦值并限制在[-1,1]范围内，避免acos产生复数
        denominator1 = 2.*BD.*L2;
        denominator2 = 2.*BD.*L23;
        
        % 避免除零错误
        denominator1(denominator1 == 0) = eps;
        denominator2(denominator2 == 0) = eps;
        
        cos_arg1 = (BD.^2+L2.^2-L4.^2)./(denominator1);
        cos_arg2 = (BD.^2+L23.^2-L3.^2)./(denominator2);
        
        % 将参数限制在有效范围内
        cos_arg1(cos_arg1 > 1) = 1;
        cos_arg1(cos_arg1 < -1) = -1;
        cos_arg2(cos_arg2 > 1) = 1;
        cos_arg2(cos_arg2 < -1) = -1;
        
        angBDC = acos(cos_arg1);
        angBDC3 = acos(cos_arg2);
        angADC = angBDC + angBDC3;
        beta = pi - angADC - theta;

        Ex = L1.*cos(theta) - L2.*cos(beta);
        
        % 确保返回实数值
        Ex(isnan(Ex)) = 0;
        Ex(isinf(Ex)) = 0;
        Ex = real(Ex);
    catch
        % 如果计算过程中出现任何错误，返回0
        Ex = zeros(size(theta));
    end
end

% 优化过程动画函数
function stop = animateOptimization(x, optimValues, state)
    stop = false;
    global optHistory iterationCount lastPlotTime
    
    switch state
        case 'iter'
            iterationCount = iterationCount + 1;
            
            % 控制绘图频率，至少间隔0.1秒才更新一次图像
            currentTime = tic;
            if toc(lastPlotTime) > 0.1
                lastPlotTime = tic;
                
                % 记录当前目标函数值
                currentFval = optimValues.fval;
                % 确保记录的是实数值
                if ~isreal(currentFval) || isnan(currentFval) || isinf(currentFval)
                    currentFval = realmax;
                end
                optHistory(end+1) = currentFval;
                
                % 更新目标函数值图
                figure(1);
                subplot(1,2,1);
                cla;
                plot(1:length(optHistory), optHistory, 'b-', 'LineWidth', 1.5);
                title('优化过程中目标函数值的变化');
                xlabel('迭代次数');
                ylabel('目标函数值');
                grid on;
                
                % 更新轨迹图
                subplot(1,2,2);
                cla;
                title(sprintf('迭代次数: %d, 当前目标函数值: %.4f', ...
                    iterationCount, currentFval));
                
                % 计算当前参数下的轨迹
                theta_plot = linspace(pi/6, pi/3, 50); % 减少点数提高动画流畅度
                params = x;
                current_trajectory = zeros(size(theta_plot));
                for i = 1:length(theta_plot)
                    current_trajectory(i) = End_Pos_calc(params(1), params(2), params(3), ...
                        params(4), params(5), theta_plot(i));
                end
                
                % 确保轨迹数据是实数
                current_trajectory = real(current_trajectory);
                current_trajectory(isnan(current_trajectory)) = 0;
                current_trajectory(isinf(current_trajectory)) = 0;
                
                plot(theta_plot, current_trajectory, 'b-', 'LineWidth', 1.5);
                hold on;
                plot(theta_plot, zeros(size(theta_plot)), 'k--', 'LineWidth', 1);
                xlabel('输入角度 \theta (弧度)');
                ylabel('末端执行器 X 坐标');
                legend('当前轨迹', '理想直线 (X=0)', 'Location', 'best');
                grid on;
                drawnow;
            end
            
        case 'interrupt'
            stop = true;
            
        case 'done'
            % 最终结果显示
            figure(1);
            subplot(1,2,2);
            title('优化完成 - 最终轨迹');
    end
end