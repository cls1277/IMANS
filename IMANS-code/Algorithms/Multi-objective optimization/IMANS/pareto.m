function [p, r] = pareto(P)
    p = [];
    r = [];
    if size(P,2)==1
        P = P';
    end
    if iscell(P)==1
        P = P{1};
    end
    N = size(P,2);
    cntn = zeros(N, 1);
    for i = 1:N
        for j = 1:N
            cnt = zeros(3,1);
            if P(i).obj(1)>P(j).obj(1)
                cnt(1) = cnt(1)+1;
            elseif P(i).obj(1)==P(j).obj(1)
                cnt(2) = cnt(2)+1;
            else
                cnt(3) = cnt(3)+1;
            end
            if P(i).obj(2)>P(j).obj(2)
                cnt(1) = cnt(1)+1;
            elseif P(i).obj(2)==P(j).obj(2)
                cnt(2) = cnt(2)+1;
            else
                cnt(3) = cnt(3)+1;
            end
            if cnt(3)==0 && cnt(2)~=2
                cntn(i) = cntn(i) + 1;
            end
        end
        pp =[P(i).obj(1);P(i).obj(2)];
        isIN = ismember(pp, p);
        if cntn(i)==0 && isIN(1)==0 && isIN(2)==0
            p = [p pp];
            r = [r P(i)];
        end
    end
end