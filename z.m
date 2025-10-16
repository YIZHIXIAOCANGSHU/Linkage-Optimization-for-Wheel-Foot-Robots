%% 四连杆末端轨迹优化
% 四连杆机构尺寸优化
clc;
clear;

% 定义机构的初始杆长
l4=0.077*sqrt(2);%初始值
l2=0.22;
l3=0.223;
l23=0.06;
l1=0.23;

initParams = [l1,l2,l3,l4,l23];

% 定义优化约束，确保杆长在合理范围内
lb = [0.16, 0.2, 0.2, 0.105, 0.05]; % 下限
ub = [0.28, 0.28, 0.28, 0.12, 0.08]; % 上限

% 定义优化选项
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');

% 执行优化
[optimizedParams, fval] = fmincon(@objectiveFun, initParams, [], [], [], [], lb, ub, [], options);

% 显示优化结果
disp('优化后的尺寸参数：');
disp(optimizedParams);
disp('优化目标函数值：');
disp(fval);

% 优化目标函数：末端轨迹尽可能接近直线
function error = objectiveFun(params)
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

    % 计算误差平方和作为目标函数
    error = sum(endEffectorX.^2);
end

% 末端位置计算函数（需根据实际四连杆结构设计）
function [Ex] = End_Pos_calc(L1, L2, L3, L4, L23, theta)
    % 实际结构需要根据四连杆几何关系进行计算
    % 四连杆结构的简单示例计算
    angBAD = pi/4 + theta;
    BD = sqrt(L2.^2+L4.^2-2.*L2.*L4.*cos(angBAD));
    angADC = acos((BD.^2+L2.^2-L4.^2)./(2.*BD.*L2))+acos((BD.^2+L23.^2-L3.^2)./(2.*BD.*L23));
    beta = pi-angADC-theta;

    Ex = L2.*cos(theta)-L1.*cos(beta);
end