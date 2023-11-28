function p = initDataSet(dataset, datalen)
    %% Get dataset name
    p = [];
    for i = 1:length(dataset)
        for j = 1:datalen(i)
            if j<10
                pth = [dataset{i} '0' num2str(j)];
            else
                pth = [dataset{i} num2str(j)];
            end
            p = [p; pth];
        end
    end
end