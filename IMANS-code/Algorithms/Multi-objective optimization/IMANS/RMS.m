function p = RMS(bkg)
    p = zeros(1, bkg.operation);
    cnt = 0;
    for i = 1:bkg.job
        for j = 1:bkg.operations(i)
            machines = cell2mat(bkg.machines{i}(j));
            cnt = cnt + 1;
            p(cnt) = machines(randi(length(machines)));
        end
    end
end