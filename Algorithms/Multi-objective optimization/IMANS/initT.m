function initT(Pop, P)
    [~, N] = size(Pop);
    Q = zeros(length(P), 1);
    for ss = 1:length(P)
        Q(ss) = sum(P(1:ss,1));
    end
    Q(length(P)) = 1;
    for i=1:N
        p = rand;
        r = find(Q>=p);
        r = r(1);
        Pop(i).add(1,2) = r;
    end
end