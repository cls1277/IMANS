function [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(p_chrom,m_chrom,FJ,bkg)
%将工序码构建成有向图的结构，有向图中包含节点集 即 所有工序
%包含边集，即两个节点之间是否含有有向边
%每个节点包含一个加工时间。每个工件的第一个工序作为开始节点的下一个节点。最多只有N个路径，N为工件个数。每个工序的最后一个工序作为结束节点的上一个工序。
%由于柔性车间调度的特殊性，导致有向图中每个节点只有两条出度两条入度。即当前工件的下一个工序，和当前机器的下一个工序。
%在这个有向图中求出，从开始节点到结束节点的最大完工时间的关键路径。
%采用的数据结构为数组加链表。数据用来存储每个工序的下标。每个数据元素后第一个元素为当前工件工序的下一个工序，第二个元素为当前机器工序的下一个工序。
%这两个元素都采用双向链表连接在数据的元素后面。实际为SH*4的矩阵。
%遍历每个工件的第一个工序，从第一个工序开始，广度遍历第二列求出一条路径，再深度遍历第三列求出一条路径。每条路径的结束标准都是，最后工序的下一列的值为0.

N = bkg.job;
H = bkg.operations';

JOBN=length(FJ);
SH=length(p_chrom);
% drawFJSP(p_chrom,m_chrom,FJ);

%每一列分布表示 当前工件工序下一个工序的索引，当前工序所在机器的下一个工序索引，当前工序工件的上一个工序的索引，当前工序所在机器的上一个工序的索引
digraph=zeros(SH,4);%由于矩阵自身的行号本身就可以表示索引值，所以省略一列
dflag=zeros(SH,4);%用于标记该点是否找到
    %先 完成解码标记每个工序的工件和工序号
    s1=p_chrom;
    s2=zeros(1,SH);
    p=zeros(1,N);
    CriticalPath=[];%关键路径，其中包含关键工序的下标
    
    for i=1:SH
        p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
        s2(i)=p(s1(i));%记录加工过程中，工件的次数
    end
    
    for i=1:SH
        t1=s1(i);%记录到当前是那个工件
        t2=s2(i);%记录当前工件是加工到第几次
        mm(i)=m_chrom(1,sum(H(1,1:t1-1))+t2);%提取该工序该次加工的机器选择，因为机器码的排列表示该工件第几次加工所选的机器，是一段表示一个工件
    end
    %构建有向图
    for i=1:SH
        to=s1(i);tm=mm(i);
        %找当前工件工序的下一个工序
        
        for j=i+1:SH
            if to==s1(j)&&dflag(i,1)==0
                digraph(i,1)=j;
                dflag(i,1)=1; %下指针已经找到
                digraph(j,3)=i;
                dflag(j,3)=1;%上指针已经找到
            end
            %找当前机器工序的下一个工序
           if tm==mm(j)&&dflag(i,2)==0
                digraph(i,2)=j;
                dflag(i,2)=1;
                digraph(j,4)=i;
                dflag(j,4)=1;
           end
            %全部都找完跳出当前循环 或者如果当前工序当前工件的是最后一个工序
           if (dflag(i,1)==1||s2(i)==H(s1(i)))&&dflag(i,2)==1
               break;
           end
           
        end
    end
    VELTime=cell(SH,2);%每个工序的最早开始时间和最晚开始时间
    e=0;
    %采用广度遍历更新每个顶点的最早开始时间和最晚开始时间,如果不采用广度优先遍历，则会在更新后续层的时候，出现前序结点来不及计算，而导致计算错误。
    %广度优先遍历的顺序是，先算所有工件的第一个工序，再是所有工件的第二个工序。
    level=zeros(1,SH);
    len=1;L=1;
    while len<=SH
        for i=1:SH
            if s2(i)==L
                level(len)=i;
                len=len+1;
            end
        end
        L=L+1;
    end
    
    %计算最早开始时间
    LastNode=e;
    for i=1:SH
       Index=level(i);
       VELTime{Index,1}=e;%VELTime{Index,2}=e; 
       t1=e;t2=e;
       if digraph(Index,3)~=0
           last=digraph(Index,3);
%            t1=time{s1(last),s2(last),mm(last)}+VELTime{last,1};
            t_ = cell2mat(bkg.times{s1(last)}(s2(last)));
            m_ = cell2mat(bkg.machines{s1(last)}(s2(last)));
            t1 = t_(m_==mm(last)) + VELTime{last,1};
       end
       
       if digraph(Index,4)~=0
           last=digraph(Index,4);
           if isempty(VELTime{last,1})             
               VELTime{last,1}=CalculateLastOperation(s1,s2,mm,digraph,last,VELTime,bkg);
           end
%            t2=time{s1(last),s2(last),mm(last)}+VELTime{last,1};
            t_ = cell2mat(bkg.times{s1(last)}(s2(last)));
            m_ = cell2mat(bkg.machines{s1(last)}(s2(last)));
            t2 = t_(m_==mm(last)) + VELTime{last,1};
       end     
       if t1>t2 %如果t1大于t2则
           VELTime{Index,1}=t1;
       else
           VELTime{Index,1}=t2;
       end
       
       if digraph(Index,1)==0 %如果是工件的最后一个工序则说明下一个节点是LastNode虚拟结束节点, 则虚拟节点的上一个节点是所有工件的最后一个工序
%            t1=time{s1(Index),s2(Index),mm(Index)}+VELTime{Index,1};%选最大的作为最后的虚拟节点
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            t1 = t_(m_==mm(Index)) + VELTime{Index,1};
           if t1>LastNode
               LastNode=t1;
           end
       end
    end
    %计算最晚开始时间
    for i=SH:-1:1
        Index=level(i); 
%         fprintf('%s %d %s %d\r\n','O',s1(Index),'.',s2(Index));
        if digraph(Index,1)==0
            if digraph(Index,2)==0
%                 VELTime{Index,2}=LastNode-time{s1(Index),s2(Index),mm(Index)}; %最后一个拓扑有序的结点的工序的最晚开始时间等于最早开始时间             
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            VELTime{Index,2}=LastNode-t_(m_==mm(Index)); %最后一个拓扑有序的结点的工序的最晚开始时间等于最早开始时间
            else
%                 t1=LastNode-time{s1(Index),s2(Index),mm(Index)};
                t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
                m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
                t1=LastNode-t_(m_==mm(Index));
                next=digraph(Index,2);
                if isempty(VELTime{next,2})
                    VELTime{next,2}=CalculateNextOperation(s1,s2,mm,digraph,next,VELTime,LastNode,bkg);
                end
%                 t2=VELTime{next,2}-time{s1(Index),s2(Index),mm(Index)};
                    t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
                    m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
                    t2 = VELTime{next,2}-t_(m_==mm(Index));
                if t1>t2 %选最小的
                    VELTime{Index,2}=t2;
                else
                    VELTime{Index,2}=t1;
                end
            end
            continue;
        end
        
        t1=e;t2=e;
        if digraph(Index,1)~=0
           next=digraph(Index,1);
%            t1=VELTime{next,2}-time{s1(Index),s2(Index),mm(Index)};%下一个节点减当前这条边的时间，但是时间在节点上所有减当前节点是时间
          t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
          m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
          t1=VELTime{next,2}-t_(m_==mm(Index));%下一个节点减当前这条边的时间，但是时间在节点上所有减当前节点是时间
          if digraph(Index,2)==0
              VELTime{Index,2}=t1;
              continue;
          end
        end
%        
        if digraph(Index,2)~=0
           next=digraph(Index,2);
           if isempty(VELTime{next,2})
               VELTime{next,2}=CalculateNextOperation(s1,s2,mm,digraph,next,VELTime,LastNode,bkg);
           end
%            t2=VELTime{next,2}-time{s1(Index),s2(Index),mm(Index)};
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            t2 = VELTime{next,2}-t_(m_==mm(Index));
        end
       if t1>t2 %选最小的
           VELTime{Index,2}=t2;
       else
           VELTime{Index,2}=t1;
       end
%        if VELTime{Index,2}<0
%            fprintf('%s %d %s %d\r\n','O',s1(Index),'.',s2(Index));
%        end
    end
    
    %用最晚开始时间减最早开始时间，剩余时间最小的结点，则为关键路径上的结点。
    Idletime=cell(SH,1);
    for i=1:SH
        Idletime{i,1}=VELTime{i,2}-VELTime{i,1};
    end
    %找到松弛时间最小的工序则为关键路径上的工序
    MinIndex=1;
    MinIdleT=Idletime{1,1};
    for i=2:SH
        if MinIdleT>Idletime{i,1}
            MinIndex=i;
            MinIdleT=Idletime{i,1};
        end
    end
    for i=1:SH
        if MinIdleT==Idletime{i,1}
            CriticalPath=[CriticalPath,i];
        end
    end
    L=length(CriticalPath);
    block=1;
    CriticalBlock(block).B=[];%表示在同一个机器上具有前后关系的顺序的 工序
    CriticalBlock(block).B=[CriticalBlock(block).B,CriticalPath(1)];m0=mm(CriticalPath(1));
    for i=2:L
        m=mm(CriticalPath(i));
        if m0==m
            CriticalBlock(block).B=[CriticalBlock(block).B,CriticalPath(i)]; %在同一个机器上相邻的工序操作为一个邻接块
        else
            block=block+1;
            CriticalBlock(block).B=[];
            CriticalBlock(block).B=[CriticalBlock(block).B,CriticalPath(i)];
        end
        m0=m;
    end
end

function [value]=CalculateLastOperation(s1,s2,mm,digraph,last,VELTime,bkg)
    Index=last;
    e=0;
    t1=e;t2=e;
       if digraph(Index,3)~=0
           zlast=digraph(Index,3);
           if isempty(VELTime{zlast,1})               
               VELTime{zlast,1}=CalculateLastOperation(s1,s2,mm,digraph,zlast,VELTime,bkg);
           end
%            t1=time{s1(zlast),s2(zlast),mm(zlast)}+VELTime{zlast,1};
            t_ = cell2mat(bkg.times{s1(zlast)}(s2(zlast)));
            m_ = cell2mat(bkg.machines{s1(zlast)}(s2(zlast)));
            t1 = t_(m_==mm(zlast)) + VELTime{zlast,1};
       end
       
       if digraph(Index,4)~=0
           zlast=digraph(Index,4);
           if isempty(VELTime{zlast,1})              
               VELTime{zlast,1}=CalculateLastOperation(s1,s2,mm,digraph,zlast,VELTime,bkg);
           end
%            t2=time{s1(zlast),s2(zlast),mm(zlast)}+VELTime{zlast,1};
            t_ = cell2mat(bkg.times{s1(zlast)}(s2(zlast)));
            m_ = cell2mat(bkg.machines{s1(zlast)}(s2(zlast)));
            t2 = t_(m_==mm(zlast)) + VELTime{zlast,1};
       end
       if t1>t2 %如果t1大于t2则
           VELTime{Index,1}=t1;
       else
           VELTime{Index,1}=t2;
       end
       value=VELTime{Index,1};
    
end

function [value]=CalculateNextOperation(s1,s2,mm,digraph,next,VELTime,LastNode,bkg)
    Index=next;
    e=0;
    t1=e;t2=e;
    if digraph(Index,1)==0
%         t1=LastNode-time{s1(Index),s2(Index),mm(Index)};
        t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
        m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
        t1 = LastNode-t_(m_==mm(Index));
    else
        znext=digraph(Index,1);
           if isempty(VELTime{znext,1})               
               VELTime{znext,1}=CalculateNextOperation(s1,s2,mm,digraph,znext,VELTime);
           end
%          t1=VELTime{znext,1}-time{s1(Index),s2(Index),mm(Index)};
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            t1 = VELTime{znext,1}-t_(m_==mm(Index));
    end
     
    if digraph(Index,2)==0
%         t2=LastNode-time{s1(Index),s2(Index),mm(Index)};
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            t2 = LastNode-t_(m_==mm(Index));
    else      
           znext=digraph(Index,2);
           if isempty(VELTime{znext,1})              
               VELTime{znext,1}=CalculateNextOperation(s1,s2,mm,digraph,znext,VELTime);
           end
%            t2=VELTime{znext,1}-time{s1(Index),s2(Index),mm(Index)};
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            t2 = VELTime{znext,1}-t_(m_==mm(Index));
     end     
       if t1>t2 %如果t1大于t2则
           VELTime{Index,1}=t2;
       else
           VELTime{Index,1}=t1;
       end
       value=VELTime{Index,1};
    
end
