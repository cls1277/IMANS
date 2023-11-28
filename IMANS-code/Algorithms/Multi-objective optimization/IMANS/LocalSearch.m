function [newDec, newAdd] = LocalSearch(offDec, PopObj, Popadd, bkg)
%     os:ind(i,bkg.job+bkg.operation+1:bkg.job+bkg.operation*2)
%     ms:ind(i,bkg.job+1:bkg.job+bkg.operation)
%     fs:ind(i,1:bkg.job)
    [N, D] = size(offDec);
    newDec = zeros(N, D);
    newAdd = zeros(N, 1);
    for i = 1:N
        o = offDec(i, bkg.job+bkg.operation+1:bkg.job+bkg.operation*2);
        m = offDec(i, bkg.job+1:bkg.job+bkg.operation);
        f = offDec(i, 1:bkg.job);
        a = [PopObj(i,:) Popadd(i,1)];
        rd = Popadd(i,2);
%         rd = getR();
        % 轮盘赌选择算子
        if rd==1
            [os, ms, fs]=LS1(o, m, f, a, bkg);
        elseif rd==2
            [os, ms, fs]=LS2(o, m, f, a, bkg);
        elseif rd==3
            [os, ms, fs]=LS3(o, m, f, bkg);
        elseif rd==4
            [os, ms, fs]=LS4(o, m, f, bkg);
        elseif rd==5
            [os, ms, fs]=LS5(o, m, f, a, bkg);
        elseif rd==6
            [os, ms, fs]=LS6(o, m, f, a, bkg);
        elseif rd==7
            [os, ms, fs]=LS7(o, m, f, a, bkg);
        end
        newDec(i,:) = [fs ms os];
        newAdd(i) = rd;
    end
    
    function r = getR()
        s = sum(local);
        Q = zeros(length(local), 1);
        for ss = 1:length(local)
            Q(ss) = sum(local(1:ss,1));
        end
        p = randi(s);
        r = find(Q>=p);
        r = r(1);
    end
end