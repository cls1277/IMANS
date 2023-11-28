function newDec = OAP(offDec, modify, bkg)
%     order after positioning
%     如果是表示MS的维度，需要判断这个machine是不是可选的，然后修改
%     如果是表示OS的维度，按照OAP的方法，但是要保证FS不变化，因为如果FS还要变化的话几乎整个决策变量都要发生变化
    [N, D] = size(offDec);
    newDec = offDec;
    dim = [modify.dim];
    %% 把能放的放进去
    for i = 1:N
        jobs = 1:bkg.job;
        cnt = bkg.operations;
        for j = 1:length(dim)
            d = dim(j);
            if d<=bkg.job
                %% FS变异的过程
                newfactory = modify(j).offdec(i);
                if newfactory<=0 || newfactory>bkg.factory
                    newfactory = randi(bkg.factory);
                end
                newDec(i,d) = newfactory;
            elseif d>bkg.job+bkg.operation
                %% OS变异的过程
                newjob = modify(j).offdec(i);
                if ismember(newjob, jobs)==0
                    newjob = jobs(randi(length(jobs)));
                end
                cnt(newjob) = cnt(newjob) - 1;
                if cnt(newjob)==0
                    jobs(jobs==newjob) = [];
                end
                newDec(i,d) = newjob;
            else
                %% MS变异的过程
                [job, idx] = getJob(d);
                newmachine = modify(j).offdec(i);
                machines = cell2mat(bkg.machines{job}(idx));
                if ismember(newmachine, machines)==0
                    newmachine = machines(randi(length(machines)));
        %           突变修复，如果当前dim的machine不在可选择内，则随机为可选择的machine
                end
                newDec(i,d) = newmachine;
            end
        end
        
        %% 把剩余的使用OAP策略放进去
        id = bkg.job+bkg.operation;
        for j = bkg.job+bkg.operation+1:bkg.job+bkg.operation*2
            if ismember(j, dim)==1
                continue;
            end
            while true
                job = offDec(i,id+1);
                if ismember(job, jobs)==1
                    id = id + 1;
                    newDec(i, j) = job;
                    cnt(job) = cnt(job) - 1;
                    if cnt(job)==0
                        jobs(jobs==job) = [];
                    end
                    break;
                end
                id = id + 1;
            end
        end
    end
    
    function [job, idx] = getJob(d)
        s = bkg.job;
        for k = 1:bkg.job
            if (s+bkg.operations(k))>=d
                job = k;
                idx = d-s;
                break;
            end
            s = s + bkg.operations(k);
        end
    end
end