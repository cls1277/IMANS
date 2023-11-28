function [p, d] = initPath(base, dataset, datalen, flag, y1, y2)
    %% Get dataset path
    p = [];
    if nargin==3
        flag = 0;
    end
    for i = 1:length(dataset)
        ds = [base dataset{i} '\'];
        if flag==0 % 遍历所有数据集
            for j = 1:datalen(i)
                if j<10
                    pth = [ds dataset{i} '0' num2str(j) '.fjs'];
                else
                    pth = [ds dataset{i} num2str(j) '.fjs'];
                end
                p = [p; pth];
            end
        elseif flag==1 % 确定使用哪个数据集
            if i==1
                j = y1;
            else
                j = y2;
            end
            if j<10
                pth = [ds dataset{i} '0' num2str(j) '.fjs'];
            else
                pth = [ds dataset{i} num2str(j) '.fjs'];
            end
            p = [p; pth];
        else % 随机一个数据集
            j = randi(datalen(i));
            if j<10
                pth = [ds dataset{i} '0' num2str(j) '.fjs'];
            else
                pth = [ds dataset{i} num2str(j) '.fjs'];
            end
            p = [p; pth];
        end
    end
    d = initDataSet(dataset, datalen);
end