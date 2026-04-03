clear; clc;
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root_dir));
addpath(fullfile(root_dir, 'src'));
addpath(fullfile(root_dir, 'benchmarks', 'cec2017'));

func_id = 10;
D = 30;
N0 = 50;
MaxFEs = 300000;

[lb, ub, ~, fobj] = Get_Functions_cec2017(func_id, D);
[BestX, BestF, HisBestFit, out] = WECO(fobj, D, N0, MaxFEs, lb, ub, struct('seed', 2026)); %#ok<ASGLU>

fprintf('CEC2017 F%d | D=%d | Best objective value: %.6e\n', func_id, D, BestF);
figure('Color','w');
plot(HisBestFit, 'LineWidth', 1.8);
xlabel('History index');
ylabel('Best-so-far fitness');
title(sprintf('WECO on CEC2017 F%d (D=%d)', func_id, D));
grid on;
