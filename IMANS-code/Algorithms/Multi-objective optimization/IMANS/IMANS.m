classdef IMANS < ALGORITHM
% <multi> <real/integer> <large/none>
% Adaptive Inverse modeling multiobjective evolutionary algorithm based on
% decomposition
% L --- 3 --- Number of mutation dimensions
% path --- 'data/bdata2.txt' --- Dataset file path
% id --- 1 --- Count of experiment
% mark --- clstcl --- Remarks of experiment
% N --- 100 --- N
% LP --- 0.5 --- LP

    methods
        function printPareto(Algorithm,Problem)
%             clc;
            fprintf('%s on %d-objective %d-variable %s (%6.2f%%), %.2fs passed...\n',...
                class(Algorithm),Problem.M,Problem.D,class(Problem),Problem.FE/Problem.maxFE*100,Algorithm.metric.runtime);
            if Problem.FE >= Problem.maxFE
                [~, p, id, mark, ~, ~] = Algorithm.ParameterSet(3,'data/bdata2.txt',1,'clstcl',100,0.5);
                fprintf('BIMMOEAD1 Done: %s %dth time end!!!\n',p(35:42),id);
%                 p = [p(1:25) 'experiment\' mark p(35:42)];
                p = [p(1:25) 'experiment\' mark '\BIMMOEAD01' p(38:42)];
                if exist(p, 'dir')==0
                    mkdir(p); % 创建文件夹
                end
                pp = p;
                p = [p '\res' num2str(id) '.txt'];
                Population = Algorithm.result{end};
                [PF, PFx] = pareto(Population);
                fout = fopen(p,'w');
                fprintf(fout,'%5.2f %6.3f\r\n',PF);
                fclose(fout);
                x = PFx.decs;
                p = [pp '\x' num2str(id) '.xlsx'];
                xlswrite(p, x, 1, 'A1');
                p = [pp '\time' num2str(id) '.txt'];
                fout = fopen(p,'w');
                fprintf(fout,'%5.2f',Algorithm.metric.runtime);
                fclose(fout);                
            end
        end
		function main(Algorithm,Problem)
            %% Experiment setting
%             Algorithm.outputFcn = @printPareto;
			[L, path, ~, ~, N, LP] = Algorithm.ParameterSet(6,'data/data.txt',1,'clstcl',80,0.6);
			%% Parameter setting
%             L = Algorithm.ParameterSet(3);
%             path = 'data/bdata.txt';
            %% Begin algorithm
            bkg = read_data(path, 2); % 2是工厂的个数
            Problem = Problem.setParms(bkg);
            Problem.N = N;
            T = ceil(Problem.N/10); % Size of neighborhood
            PP = ones(7, 1)/7;
            SM = [];
            FM = [];
%             LP = 0.5;
            
            
			%% Generate weight vectors
			[W,Problem.N] = UniformPoint(Problem.N,Problem.M);
			
			%% Detect the neighbours of each solution
			B = pdist2(W,W);
			[~,B] = sort(B,2);
			B = B(:,1:T);
			
			%% Generate random population
            randpt = randompoint(bkg, Problem.D, Problem.N);
            Population = Problem.Evaluation(randpt, zeros(Problem.N, 1),bkg);
			Z = min(Population.objs,[],1);
			
			%% Optimization 
            while Algorithm.NotTerminated(Population)
                initT(Population, PP);
%                 for i = 1:7
%                     fprintf("%f ",PP(i))
%                 end
%                 fprintf("\n");
				Offsprings=[];
				%% Modeling and reproduction
                    [~, PF] = pareto(Population);
                    idPF = ismember(Population, PF);
                    idPF0 = find(idPF==0);
                    for l = 1:length(PF)
                        if isempty(idPF0)==0
                            rd0 = randi(length(idPF0));
                            idPF(idPF0(rd0)) = 1;
                            idPF0(rd0) = [];
                        end
                    end
%                     L = 10;
                    Offspring  = Operator(Problem,Population(idPF),bkg,L);
                    % 把每个聚类都做一遍进化操作
					Offsprings = [Offsprings,Offspring];
                
                %% 执行能量调整策略
                for i = 1:length(Offsprings)
                    offsprings = Offsprings(i);
                    o = offsprings.dec(1,bkg.job+bkg.operation+1:bkg.job+bkg.operation*2);
                    m = offsprings.dec(1,bkg.job+1:bkg.job+bkg.operation);
                    f = offsprings.dec(1,1:bkg.job);
                    a = [offsprings.obj offsprings.add];
                    [os, ms, fs, as]=EnergySaveDFJSP(o, m, f, a, bkg);
                    Problem.FE = Problem.FE + 1;
                    Offsprings(i) = SOLUTION([fs ms os], as(1,1:2), 0, [as(3) offsprings.add(1,2)]);
                end
				
				%% Update the ideal point
				Z = min(Z,min(Offsprings.objs));
                ns = zeros(1, 7);
                nf = zeros(1, 7);
				for i = 1 : length(Offsprings)
					%% Global Replacement
					all_g_TCH=max(abs((Offsprings(i).obj-repmat(Z,Problem.N,1)).*W),[],2);
					best_g_TCH=min(all_g_TCH);
					Chosen_one = find(all_g_TCH(:,1)==best_g_TCH);
                    P = B(Chosen_one(1),randperm(size(B,2)));
					%% Update the solutions in P by Tchebycheff approach
					g_old = max(abs(Population(P).objs-repmat(Z,length(P),1)).*W(P,:),[],2);
					g_new = max(repmat(abs(Offsprings(i).obj-Z),length(P),1).*W(P,:),[],2);
                    better = find(g_old>=g_new,T);
                    if Offsprings(i).add(1,2)>=1 && Offsprings(i).add(1,2)<=7
                        ns(1,Offsprings(i).add(1,2)) = ns(1,Offsprings(i).add(1,2)) + length(better);
                        nf(1,Offsprings(i).add(1,2)) = nf(1,Offsprings(i).add(1,2)) + T - length(better);
                    end
					Population(P(better)) = Offsprings(i);
                end
                SM = [SM; ns];
                FM = [FM; nf];
%                 if length(SM)>=LP
                if Problem.FE>=LP*Problem.maxFE
                    for i = 1:7
                        sr = sum(SM(:,i));
                        fr = sum(FM(:,i));
                        PP(i) = sr/(sr+fr);
                    end
                    sump = sum(PP);
                    PP = PP/sump;
                    PP(PP==0) = 0.01;
                    if ~isempty(SM), SM(1, :) = [];  end
                    if ~isempty(FM), FM(1, :) = [];  end
                end
            end
		end
	end
end