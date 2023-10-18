function p = GFS(bkg)
    p = zeros(1, bkg.job);
    tf = zeros(bkg.factory, 1);
    for i = 1:bkg.job
        t = 0;
        for j = 1:bkg.operations(i)
            times = cell2mat(bkg.times{i}(j));
            t = t + min(times);
            % 假设在这个工厂内你能得到最好的资源
        end
        ids = find(tf==min(tf));
        id = randi(length(ids));
        p(1, i) = id;
        tf(id) = tf(id) + t;
    end
end