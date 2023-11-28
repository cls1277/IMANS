function GDToExcel(mark, path, cnt, alg, savePath, RP)
    savedatasets = cell(size(path,1), 1);
    for i = 1:size(path,1)
        savedatasets{i,1} = path(i,39:42);
    end
    xlswrite(savePath, savedatasets, 4, 'A2');
    xlswrite(savePath, savedatasets, 5, 'A2');
    xlswrite(savePath, alg, 4, 'B1');
    xlswrite(savePath, alg, 5, 'B1');
    
    meanEx = zeros(size(path,1), size(alg,2));
    minEx = zeros(size(path,1), size(alg,2));
    outEx = zeros(size(path,1), size(alg,2));
    for j = 1:size(path,1)
%         RP = getTruePF(cnt, path(j,:));
        ds = zeros(cnt, size(alg,2));
        PFsum = [];
        larr = zeros(2, size(path,1), cnt);
        for kk = 1:size(alg,2)
            e = alg(1, kk);
            for i = 1:cnt
                if e<10
                    p = [path(j,1:25) 'experiment\' mark '\BIMMOEAD0' num2str(e) path(j,38:42) '\res' num2str(i) '.txt'];
                else
                    p = [path(j,1:25) 'experiment\' mark '\BIMMOEAD' num2str(e) path(j,38:42) '\res' num2str(i) '.txt'];
                end
                fid = fopen(p, 'r');
                data = fscanf(fid, '%f');
                fclose(fid);
                PF = zeros(length(data)/2, 2);
                for k = 1:2:length(data)
                    PF((k+1)/2, 1) = data(k);
                    PF((k+1)/2, 2) = data(k+1);
                end
                l = size(PF, 1);
                larr(e, j,i) = l;
                PFsum = [PFsum; PF];
            end
        end
        fmax = max(PFsum,[],1)*1.1;
        fmin = max(min(PFsum,[],1),zeros(1,2));
        fmin = 0;
        [nrp,~]  = size(RP);
        ref = (RP-repmat(fmin,nrp,1))./repmat(fmax-fmin,nrp,1);
        ref = RP;
        last = 1;
        for kk = 1:size(alg,2)
            e = alg(1, kk);
            gd = zeros(cnt, 1);
            for i = 1:cnt
                now = last + larr(e, j,i) - 1;
                gd(i,1) = convergence(PFsum(last:now,:), ref, fmin, fmax);
                last = last + larr(e,j,i);
            end
            if e<10
                pa = [path(j,1:25) 'experiment\' mark '\BIMMOEAD0' num2str(e) path(j,38:42) '\gd.txt'];
            else
                pa = [path(j,1:25) 'experiment\' mark '\BIMMOEAD' num2str(e) path(j,38:42) '\gd.txt'];
            end
            hout = fopen(pa, 'w');
            fprintf(hout,'%f\t', gd');
            fprintf(hout,'\r\n');
            meanEx(j,kk) = mean(gd);
%             minEx(j,kk) = min(gd);
            fprintf(hout, "mean: %f\r\n",meanEx(j,kk));
%             fprintf(hout, "min: %f\r\n",minEx(j,kk));
            fclose(hout);
            ds(:,kk) = gd;
        end
        [~, minEx(j,:)] = sort(meanEx(j,:));
        for kk = 1:size(alg,2)
            outEx(j,minEx(j,kk)) = kk;
        end
%         p_ = ['E:\cls1277\GitHub\DGFJSP\experiment\R1\DS\GDresult' path(j,38:42) '.xlsx'];
%         xlswrite(p_, ds, 1, 'A1');
    end
    xlswrite(savePath, meanEx, 4, 'B2'); % sheet4存gd平均值
    xlswrite(savePath, outEx, 5, 'B2'); % sheet5存gd最小值
end