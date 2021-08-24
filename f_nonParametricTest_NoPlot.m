function [stats,pmask] = f_nonParametricTest_NoPlot(data1,data2)

c1 = data1; % Vector de Condicion 1
c2 = data2; %Vector de Condicion 2

data = {data1'; data2'};
paired = 'on';

method = 'bootstrap';
perms = 10000;

[stats, ~, pvals, ~] = statcond( data,'paired', paired,...
    'method',method,'naccu',perms);

[~, pmask] = fdr(pvals(:),0.05,'nonParametric');
end

