function [newp,newm,newf]=LS7(p_chrom,m_chrom,f_chrom,fitness,bkg)
    % N8
    N = bkg.job;
    SH = bkg.operation;
        s1=p_chrom;
        s2=zeros(1,SH);
        p=zeros(1,N);
        newf=f_chrom;
        for i=1:SH
            p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
            s2(i)=p(s1(i));%记录加工过程中，工件的次数
        end
    P1=[];P2=[];IP1=[];IP2=[];
    for i=1:SH
        t1=s1(i);%记录到当前是那个工件
        t2=s2(i);%记录当前工件是加工到第几次
        if f_chrom(t1)==1
            P1=[P1 p_chrom(i)]; 
            IP1=[IP1,i];
        else
            P2=[P2 p_chrom(i)];
            IP2=[IP2,i];
        end
    end
    FJ1=[];FJ2=[];
    for i=1:N
        if f_chrom(i)==1
            FJ1=[FJ1 i];
        else
            FJ2=[FJ2 i];
        end
    end

    % CriticalPath 放的是关键路径上的operator编码中的下标
    % CriticalBlock 放的是关键路径的每个块中关键步骤的下标
    % block放的是块的个数
    if fitness(1,3)==1
        [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(P1,m_chrom,FJ1,bkg); %关键路径返回的是在子染色体中关键工序的下标
    else
        [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(P2,m_chrom,FJ2,bkg);
    end
    
    % 预处理出每个机器上的工序下标
    opcnt = zeros(bkg.job, 1);
    maidx = cell(bkg.machine, 1);
    for i = 1:bkg.operation
        job = p_chrom(i);
        if f_chrom(job)~=fitness(1,3)
            continue;
        end
        opcnt(job) = opcnt(job) + 1;
        machine = m_chrom(sum(bkg.operations(1:job-1))+opcnt(job));
        maidx{machine} = [maidx{machine} i];
    end
    
    smach = 1:bkg.machine;
    for i = CriticalPath
        for m = smach
            if ismember(i, maidx{m})==1
                rd = randi(2);
                idx = find(maidx{m}==i);
                if length(maidx{m})==1
                    continue;
                end
                if (rd==1 || idx==1) && idx~=length(maidx{m}) % 设为u
                    Index1 = i;
                    Index2 = maidx{m}(ceil(rand*(length(maidx{m})-idx))+idx);
                    b = isInBlock(Index1, Index2);
                    if b~=0 % 说明都在下标是b的块里
                        BL=length(CriticalBlock(b).B);
                        if BL>1
                           % 首工序插入内部
                            Index1=1;
                            Index2=ceil(rand*(BL-2))+1;
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index1);
                                for j=Index1:Index2-1
                                    P1(j)=P1(j+1);
                                end
                                P1(Index2)=tmp;
                            else
                                tmp=P2(Index1);
                                for j=Index1:Index2-1
                                    P2(j)=P2(j+1);
                                end
                                P2(Index2)=tmp;
                            end  

                           % 尾工序插入内部
                            Index1=ceil(rand*(BL-2))+1;
                            Index2=BL; 
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index2);
                                for j=Index2:-1:Index1+1
                                    P1(j)=P1(j-1);
                                end
                                P1(Index1)=tmp;
                            else
                                tmp=P2(Index2);
                                for j=Index2:-1:Index1+1
                                    P2(j)=P2(j-1);
                                end
                                P2(Index1)=tmp;
                            end 

                          % 内部工序插入块首之前
                            Index1=1;
                            Index2=ceil(rand*(BL-2))+1;  
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index2);
                                for j=Index2:-1:Index1+1
                                    P1(j)=P1(j-1);
                                end
                                P1(Index1)=tmp;
                            else
                                tmp=P2(Index2);
                                for j=Index2:-1:Index1+1
                                    P2(j)=P2(j-1);
                                end
                                P2(Index1)=tmp;
                            end 

                          % 内部工序插入块尾之后
                            Index1=ceil(rand*(BL-2))+1;
                            Index2=BL;
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index1);
                                for j=Index1:Index2-1
                                    P1(j)=P1(j+1);
                                end
                                P1(Index2)=tmp;
                            else
                                tmp=P2(Index1);
                                for j=Index1:Index2-1
                                    P2(j)=P2(j+1);
                                end
                                P2(Index2)=tmp;
                            end      
                        end
                    else
%                         Index2 = Index2+1;
                        if fitness(1,3)==1 %update P1
                            Index1 = find(IP1==Index1);
                            Index2 = find(IP1==Index2);
                            tmp=P1(Index1);
                            for j=Index1:Index2-1
                                P1(j)=P1(j+1);
                            end
                            P1(Index2)=tmp;
                        else %update P2
                            Index1 = find(IP2==Index1);
                            Index2 = find(IP2==Index2);
                            tmp=P2(Index1);
                            for j=Index1:Index2-1
                                P2(j)=P2(j+1);
                            end
                            P2(Index2)=tmp;
                        end
                    end
                else % 设为v
                    Index2 = i;
                    Index1 = maidx{m}(ceil(rand*(idx-1)));
                    b = isInBlock(Index1, Index2);
                    if b~=0 % 说明都在下标是b的块里
                        BL=length(CriticalBlock(b).B);
                        if BL>1
                           % 首工序插入内部
                            Index1=1;
                            Index2=ceil(rand*(BL-2))+1;
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index1);
                                for j=Index1:Index2-1
                                    P1(j)=P1(j+1);
                                end
                                P1(Index2)=tmp;
                            else
                                tmp=P2(Index1);
                                for j=Index1:Index2-1
                                    P2(j)=P2(j+1);
                                end
                                P2(Index2)=tmp;
                            end  

                           % 尾工序插入内部
                            Index1=ceil(rand*(BL-2))+1;
                            Index2=BL; 
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index2);
                                for j=Index2:-1:Index1+1
                                    P1(j)=P1(j-1);
                                end
                                P1(Index1)=tmp;
                            else
                                tmp=P2(Index2);
                                for j=Index2:-1:Index1+1
                                    P2(j)=P2(j-1);
                                end
                                P2(Index1)=tmp;
                            end 

                          % 内部工序插入块首之前
                            Index1=1;
                            Index2=ceil(rand*(BL-2))+1;  
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index2);
                                for j=Index2:-1:Index1+1
                                    P1(j)=P1(j-1);
                                end
                                P1(Index1)=tmp;
                            else
                                tmp=P2(Index2);
                                for j=Index2:-1:Index1+1
                                    P2(j)=P2(j-1);
                                end
                                P2(Index1)=tmp;
                            end 

                          % 内部工序插入块尾之后
                            Index1=ceil(rand*(BL-2))+1;
                            Index2=BL;
                            Index1=CriticalBlock(b).B(Index1);
                            Index2=CriticalBlock(b).B(Index2);
                            if fitness(1,3)==1
                                tmp=P1(Index1);
                                for j=Index1:Index2-1
                                    P1(j)=P1(j+1);
                                end
                                P1(Index2)=tmp;
                            else
                                tmp=P2(Index1);
                                for j=Index1:Index2-1
                                    P2(j)=P2(j+1);
                                end
                                P2(Index2)=tmp;
                            end   
                        end                     
                    else
%                         Index1 = Index1-1;
                        if fitness(1,3)==1 %update P1
                            Index1 = find(IP1==Index1);
                            Index2 = find(IP1==Index2);
                            tmp=P1(Index2);
                            for j=Index2:-1:Index1+1
                                P1(j)=P1(j-1);
                            end
                            P1(Index1)=tmp;
                        else %update P2
                            Index1 = find(IP2==Index1);
                            Index2 = find(IP2==Index2);
                            tmp=P2(Index2);
                            for j=Index2:-1:Index1+1
                                P2(j)=P2(j-1);
                            end
                            P2(Index1)=tmp;
                        end
                    end
                end 
                smach(smach==m) = [];
                break;
            end
        end
    end
    
newm=m_chrom;

newp=zeros(1,SH);
L=length(IP1);
for i=1:L
  newp(1,IP1(i))=P1(i);  
end

L=length(IP2);
for i=1:L
  newp(1,IP2(i))=P2(i);  
end
    
    function p = isInBlock(idx1, idx2)
        p = 0;
        for k = 1:block
            B = CriticalBlock(k).B;
            if ismember(idx1,B)==1 && ismember(idx2,B)==1
                p = k;
                break;
            end
        end
    end
end