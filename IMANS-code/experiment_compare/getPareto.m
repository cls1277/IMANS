%% Get the set of pareto
function getPareto(mark, cnt, path, alg)
    for j = 1:size(path,1)
        ds = zeros(cnt, size(alg,2));
        PFsum = [];
        larr = zeros(2, size(path,1), cnt);
        for e = alg
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
                PF = unique(PF, 'rows');
                l = size(PF, 1);
                larr(e, j,i) = l;
                PFsum = [PFsum; PF];
            end
        end
        %% Calculate HV
        fmin = max(min(PFsum,[],1),zeros(1,2));
        fmax = max(PFsum,[],1)*1.1;
        last = 1;   
        ec = 0;
        for e = alg
            ec = ec + 1;
            hvsum = zeros(1, cnt);
            for i = 1:cnt
                now = last + larr(e, j,i) - 1;
                hv = getHV(PFsum(last:now,:), fmin, fmax);
                hvsum(1,i) = hv;
                last = last + larr(e, j,i);
            end
            if e<10
                p = [path(j,1:25) 'experiment\' mark '\BIMMOEAD0' num2str(e) path(j,38:42) '\hv.txt'];
            else
                p = [path(j,1:25) 'experiment\' mark '\BIMMOEAD' num2str(e) path(j,38:42) '\hv.txt'];
            end
            hout = fopen(p, 'w');
            fprintf(hout,'%f\t', hvsum);
            fprintf(hout,'\r\n');
            fprintf(hout, "mean: %f\r\n",mean(hvsum));
            fprintf(hout, "max: %f\r\n",max(hvsum));
            fclose(hout);
            ds(:,ec) = hvsum';
        end
%         p_ = ['E:\cls1277\GitHub\DGFJSP\experiment\R1\DS\HVresult' path(j,38:42) '.xlsx'];
%         xlswrite(p_, ds, 1, 'A1');
    end
end