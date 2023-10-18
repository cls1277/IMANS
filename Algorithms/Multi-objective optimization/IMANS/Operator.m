function [Offspring, local] = Operator(Problem,Population,bkg,L)
% The Gaussian process based reproduction

%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is modified from the code in
% http://www.soft-computing.de/jin-pub_year.html

    %% Parameter setting
    if nargin < 3
        L = 3;
    end
    PopDec = Population.decs;
	PopObj = Population.objs;
    PopAdd = Population.adds;
    [N,D]  = size(PopDec);
    
    %% 使用局部搜索算子
    OffDec = [];
    OffAdd = [];
    [offDec, offAdd] = LocalSearch(PopDec, PopObj, PopAdd, bkg);
    OffDec = [OffDec;offDec];
    OffAdd = [OffAdd;offAdd];
%     for i = 1:length(local)
%         fprintf("%d ",local(i));
%     end
%     fprintf("\n");
    
    %% Gaussian process based reproduction
    if length(Population) < 2*Problem.M
        OffDec = PopDec;
    else
        fmin   = 1.5*min(PopObj,[],1) - 0.5*max(PopObj,[],1);
        fmax   = 1.5*max(PopObj,[],1) - 0.5*min(PopObj,[],1);
        % Train one groups of GP models for each objective
        for m = 1 : Problem.M
            parents = randperm(N,floor(N/Problem.M)); % 在这个subpopulation中选出这几个个体
%             parents = randperm(N, N);
            offDec = PopDec(parents,:); % 取出这几个个体的决策变量
            MODY = struct('dim', [], 'offdec', []);
            modify = repmat(MODY, L, 1);
            cnt = 0;
            for d = randperm(D,L) % 突变决策变量的第d维上的数据
                % Gaussian Process
                cnt = cnt + 1;
                modify(cnt).dim = d;
                try
                    [ymu,ys2] = gp(struct('mean',[],'cov',[],'lik',log(0.01)),...
                                   @infExact,@meanZero,@covLIN,@likGauss,...
                                   PopObj(parents,m),PopDec(parents,d),...
                                   linspace(fmin(m),fmax(m),size(offDec,1))');
                      modify(cnt).offdec = round(ymu + rand*sqrt(ys2).*randn(size(ys2)));
                catch
                end
            end
            offDec = AfterGP(offDec, modify, bkg);
            OffDec = [OffDec;offDec];
            OffAdd = [OffAdd;zeros(size(offDec,1),1)];
        end
    end
   
    Offspring = Problem.Evaluation(OffDec, OffAdd, bkg);
end