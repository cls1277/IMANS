function p = GMS2(bkg, fs, os)
    % 选那个最小能耗的机器
    E = [4; 1];
    time_m = zeros(bkg.factory, bkg.machine);
    time_p = zeros(bkg.factory, bkg.machine); % 每台机器的运行时间
    time_i = zeros(bkg.factory, bkg.machine); % 每台机器的空闲时间
    time_j = zeros(bkg.job, 1);
    p = zeros(1, bkg.operation);
    cnt = zeros(bkg.job, 1);
    for i = 1:length(os)
        job = os(i);
        factory = fs(job);
        cnt(job) = cnt(job) + 1;
        machines = cell2mat(bkg.machines{job}(cnt(job)));
        times = cell2mat(bkg.times{job}(cnt(job)));
        p_se = zeros(length(machines), 1);
        for j = 1:length(machines)
            machine = machines(j);
            time = times(j);
            t = max(time_m(factory, machine), time_j(job));
            t_p = time_p(factory, machine) + time;
            t_i = time_i(factory, machine) + time_m(factory, machine) - t;
            p_se(j) = t_p*E(1) + t_i*E(2);
        end
        [~, idx] = min(p_se);
        p(1,op2dim(job,cnt(job))) = machines(idx);
        time = times(idx);
        time_j(job) = time_j(job) + time;
        time_m(factory, machine) = time_m(factory, machine) + time;
    end    
    
    function y = op2dim(j, c)
        y = sum(bkg.operations(1:j-1)) + c;
    end
end