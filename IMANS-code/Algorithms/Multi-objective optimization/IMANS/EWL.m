function newDec = EWL(offDec, modify, bkg)
%     exchange with locations
%     又分三种情况，第1个location 最后一个location 和 随机一个location
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
                ids = find(newDec(i,bkg.job+bkg.operation+1:D)==newjob);
                rd = randi(3);
                if rd==1
                    id = ids(1);
                elseif rd==2
                    id = ids(length(ids));
                else
                    id = ids(randi(length(ids)));
                end
                id = id + bkg.job + bkg.operation;
                [newDec(i,id), newDec(i,d)] = swap(newDec(i,id), newDec(i,d));
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
    end
    
    function [a, b] = swap(a, b)
        c = a;
        a = b;
        b = c;
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