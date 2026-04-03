clear; clc;
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(root_dir));
addpath(fullfile(root_dir, 'src'));

D = 30;
N0 = 50;
MaxFEs = 10000;
lb = -100;
ub = 100;

obj_fun = @(x) sum(x.^2);
[BestX, BestF, HisBestFit, out] = WECO(obj_fun, D, N0, MaxFEs, lb, ub, struct('seed', 2026, 'verbose', true)); %#ok<ASGLU>

fprintf('Best objective value: %.6e\n', BestF);
figure('Color','w');
plot(HisBestFit, 'LineWidth', 1.8);
xlabel('History index');
ylabel('Best-so-far fitness');
title('WECO on Sphere');
grid on;
