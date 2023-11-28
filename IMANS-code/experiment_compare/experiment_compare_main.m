function experiment_compare_main()
    %% Dataset setting
    dataset = {'Mk', 'DP'};
    datalen = [10, 10];
%     dataset = {'RR'};
%     datalen = [1];
    base = 'C:\cls\1\cls1277\GitHub\DGFJSP\instances\';
    
    %% Experiment setting
    [paths, datasets] = initPath(base, dataset, datalen); % 数据集路径
    mark = 'return'; % 实验备注
    cnt = 20; % 实验次数
    iter = 65000; % 最大评价次数
    path = paths;

    %% Do experiment
    Lset = [3 6 9 11];
    Nset = [40 60 80 100];
    LPset = [0.4 0.5 0.6 0.7];
    table = [2,3,3]; % 通过设置table，可以设置不同的参数组合
    for i = 1:1
        L = Lset(1,table(i,1));
        N = Nset(1,table(i,2));
        LP = LPset(1,table(i,3));
        experiment_compare_1(path, mark, cnt, iter, L, N, LP);
    end
    return ;
    %% 
    % 下面是hv gd spread相关的，路径是绝对路径，
    % 是当时比较评价指标时候的代码，仅供参考
    % 在不修改路径的情况下可以直接看到如何处理的实验结果和用指标评估的
        
    %% Hv to excel(1 2)
    getPareto(mark, cnt, path, alg);
    pathHV = initHVPath(base1, dataset, datalen, mark, alg);
    hvToExcel(pathHV, resultPath, datasets, alg);
    fprintf("hv to excel done...\n");
%     return ;

    %% GD to excel(4 5)
    ref_point=[[0,0];[5,0];[10,0];[0,5];[0,10]];
    GDToExcel(mark, path, cnt, alg, resultPath, ref_point);
    fprintf("gd to excel done...\n");
%     return ;

    %% Spread to excel(6 7)
    SpreadToExcel(mark, path, cnt, alg, resultPath, ref_point);
    fprintf("spread to excel done...\n");
%     return ; 
end