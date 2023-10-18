function p = GMS1(bkg, fs, os)
    % 选那个最小时间的机器
    time_m = zeros(bkg.factory, bkg.machine);
    time_j = zeros(bkg.job, 1);
    p = zeros(1, bkg.operation);
    cnt = zeros(bkg.job, 1);
    for i = 1:length(os)
        job = os(i);
        factory = fs(job);
        cnt(job) = cnt(job) + 1;
        machines = cell2mat(bkg.machines{job}(cnt(job)));
        times = cell2mat(bkg.times{job}(cnt(job)));
        t_se = zeros(length(machines), 1);
        for j = 1:length(machines)
            machine = machines(j);
            t_se(j) = max(time_m(factory, machine), time_j(job));
        end
        [~, idx] = min(t_se);
        p(1,op2dim(job,cnt(job))) = machines(idx);
        time = times(idx);
        time_j(job) = time_j(job) + time;
        time_m(factory, machine) = time_m(factory, machine) + time;
    end
    
    function y = op2dim(j, c)
        y = sum(bkg.operations(1:j-1)) + c;
    end
end