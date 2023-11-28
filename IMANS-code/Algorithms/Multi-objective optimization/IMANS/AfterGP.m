function newDec = AfterGP(offDec, modify, bkg)
    % 有两种处理方式：除了modify内的东西都按顺序放；和一个位置互换
    rd = randi(2);
    if rd == 1
        newDec = OAP(offDec, modify, bkg);
    else
        newDec = EWL(offDec, modify, bkg);
    end
end