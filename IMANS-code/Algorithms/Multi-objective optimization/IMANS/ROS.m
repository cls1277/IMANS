function p = ROS(bkg)
    cnt = bkg.operations;
    jobs = 1:bkg.job;
    p = zeros(1, bkg.operation);
    for i = 1:bkg.operation
        job = jobs(randi(length(jobs)));
        p(i) = job;
        cnt(job) = cnt(job) - 1;
        if cnt(job)==0
            jobs(jobs==job) = [];
        end
    end
end