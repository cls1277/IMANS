function ind = randompoint(bkg, d, n)
    %RANDOMNEW to generate n new point randomly from the mop problem given.
    ind = zeros(n, d);
    % 3行1列的矩阵用于随机决定使用哪个初始化算子
    for i = 1:n
        %% 初始化FS
        r = round(rand*1);
        if r==0
            fs = RFS(bkg);
        else
            fs = GFS(bkg);
        end
        ind(i,1:bkg.job) = fs;
        %% 初始化OS
        os = ROS(bkg);
        ind(i,bkg.job+bkg.operation+1:bkg.job+bkg.operation*2) = os;
        %% 初始化MS
        r = round(rand*2);
        if r==0
            ms = RMS(bkg);
        elseif r==1
            ms = GMS1(bkg, fs, os);
        else
            ms = GMS2(bkg, fs, os);
        end
        ind(i,bkg.job+1:bkg.job+bkg.operation) = ms;
    end
end 
