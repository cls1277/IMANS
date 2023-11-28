function bkg = read_data(path, factorys)
    fid = fopen(path, 'r');
    data = fscanf(fid, '%d');
    bkg = struct('factory', [], 'job', [], 'machine', [], 'operation', [], 'operations', [], 'machines', [], 'times', []);
    bkg.factory = factorys;
    bkg.job = data(1);
    bkg.machine = data(2);
    % data(3)没用到
    bkg.operations = zeros(bkg.job, 1);
    bkg.machines = cell(bkg.job, 1);
    bkg.times = cell(bkg.job, 1);
    cnt = 3;
    for i = 1:bkg.job
        cnt = cnt + 1;
        bkg.operations(i) = data(cnt);
        operations_m = cell(bkg.operations(i), 1);
        operations_t = cell(bkg.operations(i), 1);
        for j = 1:bkg.operations(i)
            cnt = cnt + 1;
            num = data(cnt);
            machines = zeros(1, num);
            times = zeros(1, num);
            for k = 1:num
                cnt = cnt + 1;
                machine = data(cnt);
                cnt = cnt + 1;
                time = data(cnt);
                machines(k) = machine;
                times(k) = time;
            end
            operations_m{j} = machines;
            operations_t{j} = times;
        end
        bkg.machines{i} = [bkg.machines{i} operations_m];
        bkg.times{i} = [bkg.times{i} operations_t];
    end
    bkg.operation = sum(bkg.operations);
end