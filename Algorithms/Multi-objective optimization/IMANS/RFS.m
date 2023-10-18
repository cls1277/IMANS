function p = RFS(bkg)
    p = round(rand(1,bkg.job)*(bkg.factory-1)+1);
end