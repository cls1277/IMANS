% if count(py.sys.path,'') == 0
%     insert(py.sys.path, int32(0), '');
% end
% 
% % fid = fopen('moon.txt','r');
% % data = fscanf(fid, '%f');
% % x = py.list({});
% % for i = 1:2:length(data)
% %     l = py.list({data(i), data(i+1)});
% %     x.append(l);
% % end
% x = rand(100, 2);
% % py.print(x);
% x = py.numpy.array(x);
% pyModule = py.importlib.import_module('FDP');
% py.importlib.reload(pyModule);
% fdp = pyModule.FastDensityPeaks();
% labels = fdp.fit_predict(x);
% % labels = double(labels.tolist());
% ls = py.array.array('d',labels.tolist());
% p = double(ls);
% py.print(labels);
% % fprintf("hello world\n");
% 
% a = rand(5, 1)
% a(3) = [];
% a(3)
% a(4)
% 
% fid = fopen('test.txt','r');
% data = fscanf(fid, '%f');
% x = zeros(68, 1);
% y = zeros(68, 1);
% for i = 1:2:length(data)
%     x((i+1)/2) = data(i);
%     y((i+1)/2) = data(i+1);
% end
% [ymu,ys2] = gp(struct('mean',[],'cov',[],'lik',log(0.01)), @infExact,@meanZero,@covLIN,@likGauss,x,y,linspace(482448,643596,34)');
% 
% round(rand(1,4)*(2-1)+1)

% a = rand(2,3)
% max(a(:))
% 
% x = gpml_randn(0.8, 20, 1);                 % 20 training inputs
% y = sin(3*x) + 0.1*gpml_randn(0.9, 20, 1);  % 20 noisy training targets
% xs = linspace(-3, 3, 61)';                  % 61 test inputs 
% meanfunc = [];                    % empty: don't use a mean function
% covfunc = @covSEiso;              % Squared Exponental covariance function
% likfunc = @likGauss;              % Gaussian likelihood
% hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
% hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
% [mu,s2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, x, y, xs);
% f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
% fill([xs; flipdim(xs,1)], f, [7 7 7]/8);
% hold on; plot(xs, mu); plot(x, y, '+');

% a = [1, 2, 2, 2, 3];
% find(a==2)

% round(rand*1)
% a = rand(2, 3)
% [~, b] = max(a(:))

% local = [4,2,3,4]';
% s = sum(local);
% Q = zeros(length(local), 1);
% for ss = 1:length(local)
%     Q(ss) = sum(local(1:ss,1));
% end
% while true
%     p = randi(s);
%     r = find(Q>=p);
%     r = r(1);
%     fprintf("%d\n", r);
%     if r==0
%         break;
%     end
% end
% 
% p = 'E:\cls1277\GitHub\DGFJSP\instances\DP\DP01.fjs';
% % dataset = p(35:42)
% p = [p(1:25) 'experiment' p(35:42) '\res' num2str(2) '.txt']

% p = 'E:\cls1277\GitHub\DGFJSP\test\test.txt';
% f = fopen(p, 'w');
% h = rand(5, 1);
% fprintf(f, "%f\r\n", h);
% fprintf(f, "mean: %f\r\n",mean(h));
% fprintf(f, "max: %f\r\n",max(h));
% fclose(f);

% rand

numst = 7;
ps = 100;
pfit = ones(1, numst);
normfit=zeros(1, numst);
rr = rand;
spacing = 1/ps;
randnums = sort(mod(rr : spacing : 1 + rr - 0.5 * spacing, 1));

sumfit=0;
for j=1:numst%算子概率归一化
  sumfit=pfit(j)+sumfit;
end
for j=1:numst
  normfit(j)=pfit(j)/sumfit;
end
partsum = 0;
count(1) = 0;
stpool = [];

for j = 1 : length(pfit)
 partsum = partsum + normfit(j);
 count(j + 1) = length(find(randnums < partsum));
 select(j, 1) = count(j + 1) - count(j);
 stpool = [stpool; ones(select(j, 1), 1) * j];
end
stpool = stpool(randperm(ps));
size(stpool)