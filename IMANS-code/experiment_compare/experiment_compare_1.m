function experiment_compare_1(path, mark, cnt, iter, L, N, LP)
    %% Begin experiment
    for j = 1:size(path,1)
        for i = 1:cnt
            algSet = cell(1,3);
            algSet{1,1} = @IMANS;
            algSet{1,2} = L;
            algSet{1,3} = path(j,:);
            algSet{1,4} = i;
            algSet{1,5} = mark;
            algSet{1,6} = N;
            algSet{1,7} = LP;
            %% do it
            platemo('algorithm',algSet,'problem',@DGFJSP, 'maxFE', iter);
        end
    end
    fprintf("experiment 1 end...\n");
end

