classdef DGFJSP < PROBLEM
% <multi> <integer> <large/none> <expensive/none>
    methods
        %% Default settings of the problem(false)
        function Setting(obj)
            obj.M = 2;
        end
        %% Initialization of vectors(false)
        function Population = Initialization(N)
        end
        %% Calculate objective values
        function Population = Evaluation(obj, randpt, randadd, bkg)
            [N, ~] = size(randpt);
            Population = repmat(SOLUTION, 1, N);
            for i = 1:N
                z = evaluate(randpt(i, :), bkg);
                % 取log是防止值太大
%                 Population(1, i) = SOLUTION(randpt(i, :), [z(1) log(z(2))], 0, z(3));
                Population(1, i) = SOLUTION(randpt(i, :), [z(1) z(2)], 0, [z(3) randadd(i)]);
                % add中存的值是关键工厂的下标
                % add的第二个值表示要选择哪个算子
                % add的第三个值表示这个决策变量离哪个权重向量近（已证明不优）
            end
            obj.FE = obj.FE + N;
        end
        %% Default settings of the problem(true)
        function p = setParms(obj, bkg)
            p = obj;
            p.M = 2;
            p.D = bkg.job+2*bkg.operation;
            p.encoding = ones(1, p.D)*2;
            p.lower = [ones(1,bkg.job) ones(1,bkg.operation) ones(1,bkg.operation)];
            p.upper = [ones(1,bkg.job)*bkg.factory ones(1,bkg.operation)*bkg.machine ones(1,bkg.operation)*bkg.job];
        end
    end
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