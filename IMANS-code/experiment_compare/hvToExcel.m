function hvToExcel(paths, savePath, datasets, algs)
%     datasets是使用的数据集的名字的列向量
%     algs是使用的算法的名字的列向量
    savedatasets = cell(size(datasets,1), 1);
    for i = 1:size(datasets,1)
        savedatasets{i,1} = datasets(i,:);
    end
    xlswrite(savePath, savedatasets, 1, 'A2');
    xlswrite(savePath, savedatasets, 2, 'A2');
    xlswrite(savePath, algs, 1, 'B1');
    xlswrite(savePath, algs, 2, 'B1');

    meanEx = zeros(size(datasets,1), size(algs,1));
    maxEx = zeros(size(datasets,1), size(algs,1));
    for i = 1:size(paths,1)
        path = paths(i,:);
        [dataset, alg] = solvePath(path);
        for k = 1:size(datasets, 1)
            if strcmp(datasets(k,:), dataset)
                dataset = k;
                break;
            end
        end
        alg = find(algs==str2double(alg));
        fid = fopen(path);
        lines = fgetl(fid);
        while ischar(lines)
            if contains(lines,'mean')
                meanHV = str2double(lines(1, 6:size(lines,2)));
            end
            if contains(lines,'max')
                maxHV = str2double(lines(1, 5:size(lines,2)));
            end            
            lines = fgetl(fid);
        end
        fclose(fid);
        meanEx(dataset, alg) = meanHV;
        maxEx(dataset, alg) = maxHV;
    end
    outEx = zeros(size(maxEx));
    for k = 1:size(datasets, 1)
        [~, maxEx(k,:)] = sort(meanEx(k,:), 'descend');
        for kk = 1:size(maxEx,2)
            outEx(k,maxEx(k,kk)) = kk;
        end
    end
    xlswrite(savePath, meanEx, 1, 'B2'); % sheet1存hv平均值
    xlswrite(savePath, outEx, 2, 'B2'); % sheet2存hv最大值

    function [d, a] = solvePath(path)
        d = [];
        a = [];
        idx = strfind(path, 'Mk');
        if isempty(idx)
            idx = strfind(path, 'DP');
        end
        for j = idx:size(path,2)
            if strcmp(path(1,j), '\')==1
                break;
            end
            d = [d path(1,j)];
        end
        idx = strfind(path, 'BIMMOEAD') + 8;
        for j = idx:size(path,2)
            if strcmp(path(1,j), '\')==1
                break;
            end
            a = [a path(1,j)];
        end
    end
end