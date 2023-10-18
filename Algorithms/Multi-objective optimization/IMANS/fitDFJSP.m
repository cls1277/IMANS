function [fit1,fit2,fit3]=fitDFJSP(p_chrom,m_chrom,f_chrom,bkg)
    z = evaluate([f_chrom m_chrom p_chrom], bkg);
    fit1 = z(1);
    fit2 = z(2);
    fit3 = z(3);
end

function p = evaluate(x, bkg)
    % 常数：Ep加工功率(4) Es等待功率(1)
    E = [4; 1];
    if size(x,2)==1
        x = x';
    end
    % x是行向量
    p = zeros(3, 1);
    time_m = zeros(bkg.factory, bkg.machine);
    time_j = zeros(bkg.job, 1);
    cnt_j = zeros(bkg.job, 1);
    time_E = zeros(1, 2);
    for i = bkg.job+bkg.operation+1:bkg.job+bkg.operation*2
        job = x(i);
        factory = x(job);
        cnt_j(job) = cnt_j(job) + 1;
        machine = getMachine(job, cnt_j(job));
        machines = cell2mat(bkg.machines{job}(cnt_j(job)));
        times = cell2mat(bkg.times{job}(cnt_j(job)));
        time = times(machines==machine);
        t = max(time_m(factory, machine), time_j(job));
        time_E(1) = time_E(1) + time;
        time_E(2) = time_E(2) + t - time_m(factory, machine);
        time_m(factory, machine) = t+time;
        time_j(job) = t+time;
    end
    [p(1), idx] = max(time_m(:));
    p(2) = time_E*E;
    p(3) = mod(idx, bkg.factory);
    if p(3)==0
        p(3) = bkg.factory;
    end
    
    function m = getMachine(j, c)
        pre = bkg.job + sum(bkg.operations(1:j-1)) + c;
        m = x(pre);
    end
end