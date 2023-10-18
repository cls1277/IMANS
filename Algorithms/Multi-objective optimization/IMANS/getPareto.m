path = 'Data/AIMMOEAD/AIMMOEAD_TA01_M2_D32_32.mat';
load(path);
solution = result(size(result,1),2);
p = pareto(solution);
s = [p.obj];
x = s(1:2:length(s));
y = s(2:2:length(s));
scatter(x, y);