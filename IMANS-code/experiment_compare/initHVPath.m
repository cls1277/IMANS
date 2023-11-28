function p = initHVPath(base, dataset, datalen, mark, alg)
    %% Get hv file path
    p = [];
    for k = alg
        for i = 1:length(dataset)
            if k<10
                ds = [base mark '\BIMMOEAD0' num2str(k) '\'];
            else
                ds = [base mark '\BIMMOEAD' num2str(k) '\'];
            end
            for j = 1:datalen(i)
                if j<10
                    pth = [ds dataset{i} '0' num2str(j) '\hv.txt'];
                else
                    pth = [ds dataset{i} num2str(j) '\hv.txt'];
                end
                p = [p; pth];
            end
        end
    end
end