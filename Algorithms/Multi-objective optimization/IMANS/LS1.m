function [newp,newm,newf]=LS1(p_chrom,m_chrom,f_chrom,fitness,bkg)
    % N5
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
    for i = 2:block-1
        BL=length(CriticalBlock(i).B);
        if BL>1
            Index1=BL-1;
            Index2=BL;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            if fitness(1,3)==1
                tmp=P1(Index1);
                P1(Index1)=P1(Index2);
                P1(Index2)=tmp;
            else
                tmp=P2(Index1);
                P2(Index1)=P2(Index2);
                P2(Index2)=tmp;
            end
            
            Index1=1;
            Index2=2;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            if fitness(1,3)==1
                tmp=P1(Index1);
                P1(Index1)=P1(Index2);
                P1(Index2)=tmp;
            else
                tmp=P2(Index1);
                P2(Index1)=P2(Index2);
                P2(Index2)=tmp;
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

end