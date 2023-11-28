function [p_chrom,m_chrom,f_chrom,fitness]=EnergySaveDFJSP(p_chrom,m_chrom,f_chrom,fitness,bkg) %通过减少等待时间来减少能量消耗
%     global N SH AP AM AN AF;
    N = bkg.job;
    SH = bkg.operation;

    s1=p_chrom;
    s2=zeros(1,SH);
    p=zeros(1,N);
    
    for i=1:SH
        p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
        s2(i)=p(s1(i));%记录加工过程中，工件的次数
    end
    
P1=[];P2=[];IP1=[];IP2=[];
for i=1:SH
    t1=s1(i);%记录到当前是那个工件
    t2=s2(i);%记录当前工件是加工到第几次
    if f_chrom(t1)==1
        % P1存0工厂的所有工序
        % IP1存0工厂的所有工序的下标
        P1=[P1 p_chrom(i)];
        IP1=[IP1,i];
    else
        P2=[P2 p_chrom(i)];
        IP2=[IP2,i];
    end
end
% FJ存该工厂有哪几个工件
FJ1=[];FJ2=[];
for i=1:N
    if f_chrom(i)==1
        FJ1=[FJ1 i];
    else
        FJ2=[FJ2 i];
    end
end

[P1]=SAS2AS(P1,m_chrom,FJ1,bkg);%主动解码
[P2]=SAS2AS(P2,m_chrom,FJ2,bkg);

[P1]=AS2FAS(P1,m_chrom,FJ1,bkg);%全主动解码
[P2]=AS2FAS(P2,m_chrom,FJ2,bkg);
new_f=f_chrom;
new_m=m_chrom;

% 把解码后的工厂的工序按照下标放到一个新的向量里
% 此时工厂码和机器码不发生变化
new_p=zeros(1,SH);
L=length(IP1);
for i=1:L
  new_p(1,IP1(i))=P1(i);  
end

L=length(IP2);
for i=1:L
  new_p(1,IP2(i))=P2(i);  
end

% new_p2=zeros(1,SH);
% L=length(IP1);
% for i=1:L
%   new_p2(1,IP1(i))=P3(i);  
% end
% 
% L=length(IP2);
% for i=1:L
%   new_p2(1,IP2(i))=P4(i);  
% end

[newfit(1,1),newfit(1,2),newfit(1,3)]=fitDFJSP(new_p,new_m,new_f,bkg);
% [newfit2(1,1),newfit2(1,2),newfit2(1,3)]=fitDFJSP(new_p2,new_m,new_f);        
            if NDS(newfit,fitness)==1 %新解支配旧解
                p_chrom=new_p;
                m_chrom=new_m;
                f_chrom=new_f;
                fitness=newfit;           
            end
end

function [newp]=SAS2AS(p_chrom,m_chrom,FJ,bkg)
%     global N H TM time;
    N = bkg.job;
    H = bkg.operations';
    TM = bkg.machine;
    e=0;
    JOBN=length(FJ);
    SH=length(p_chrom);
    finish={};%工序完工时间
    start={};%工序开始时间
    for i=1:JOBN%初始化完成时间矩阵
        JOBI=FJ(i);
        for j=1:H(JOBI)
            finish{JOBI,j}=e;
            start{JOBI,j}=e;
        end
    end
    
    mt=cell(1,TM);
    for i=1:TM%初始化机器最大完成时间数组
        Machine(i).Op=[];
        Machine(i).GapT=[];
        Machine(i).MFT=[];
        mt{i}=e;
    end
    s1=p_chrom;
    s2=zeros(1,SH);
    p=zeros(1,N);
    
    for i=1:SH
        p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
        s2(i)=p(s1(i));%记录加工过程中，工件的次数
    end
    
    for i=1:SH
        t1=s1(i);%记录到当前是那个工件
        t2=s2(i);%记录当前工件是加工到第几次
        mm(i)=m_chrom(1,sum(H(1,1:t1-1))+t2);
        %提取该工序该次加工的机器选择，因为机器码的排列表示该工件第几次加工所选的机器，是一段表示一个工件
    end
    %开始解码
    for i=1:SH
        if(s2(i)==1)
         ON=length(Machine(mm(i)).Op);%该机器上目前的工序数   
         if ON>0
%             t=time{s1(i),s2(i),mm(i)};
            t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
            m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
            t = t_(m_==mm(i));
            Index1=0;
            for j=1:ON
                if Machine(mm(i)).GapT(j)-t>0
                    Index1=j; %在Index工序前插入新来的工序
                    break;
                end
            end
            if Index1~=0
                Index1=Machine(mm(i)).Op(Index1);% 在s1向量中把s1(i)插入到s1(Index1)前面,Index1<i
                tmp=s1(i);
                for j=i:-1:Index1+1
                    s1(j)=s1(j-1);
                end
                s1(Index1)=tmp;
                tmp=s2(i);
                for j=i:-1:Index1+1
                    s2(j)=s2(j-1);
                end
                s2(Index1)=tmp;
                
                tmp=mm(i);
                for j=i:-1:Index1+1
                    mm(j)=mm(j-1);
                end
                mm(Index1)=tmp;
                
            for j=1:ON
                if Machine(mm(Index1)).Op(j)>=Index1
                    Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)+1;
                    %由于前面插入工序，其他索引需要整体向后退一格
                end
            end
            for k=1:TM
                if k~=mm(Index1)
                    ON2=length(Machine(k).Op);
                    for h=1:ON2
                        if Machine(k).Op(h)>Index1&&Machine(k).Op(h)<i
                            Machine(k).Op(h)=Machine(k).Op(h)+1;
                        end
                    end
                end
            end
            Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
            tmp2=Machine(mm(Index1)).Op;
            Machine(mm(Index1)).Op=sort(tmp2,'ascend');
            IIndex=find(Machine(mm(Index1)).Op==Index1);
            if IIndex==1
                start{s1(Index1),s2(Index1)}=0;
            else
                LastOp=Machine(mm(Index1)).Op(IIndex-1);
                start{s1(Index1),s2(Index1)}=max(0,finish{s1(LastOp),s2(LastOp)});
            end   
                finish{s1(Index1),s2(Index1)}=t+start{s1(Index1),s2(Index1)};
                ON=ON+1;
            for j=1:ON
                Index1=Machine(mm(Index1)).Op(j);
                if j==1
                    Machine(mm(Index1)).GapT(j)=0;                   
                else
                    LastOp=Machine(mm(Index1)).Op(j-1);
                    Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                    if Machine(mm(Index1)).GapT(j)<0
                        Machine(mm(Index1)).GapT(j)
                    end
                end
                Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
            end
            mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else %如果index1==0说明没有空位就要老老实实的去在后面
                start{s1(i),s2(i)}=Machine(mm(i)).MFT(ON);
%                 mt{mm(i)}=start{s1(i),s2(i)}+time{s1(i),s2(i),mm(i)}; 
                t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
                m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
                mt{mm(i)} = start{s1(i),s2(i)}+t_(m_==mm(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
            
         else
%              mt{mm(i)}=time{s1(i),s2(i),mm(i)};
             t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
             m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
             mt{mm(i)} = t_(m_==mm(i));
             start{s1(i),s2(i)}=0;
             finish{s1(i),s2(i)}=mt{mm(i)};
             Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
             Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];
             Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
         end
         
        else
            %开始在相同的机器上寻找合适的空位插入，如果自身的加工时间小于空闲时间则可以插入
         ON=length(Machine(mm(i)).Op);%该机器上目前的工序数   
         if ON>0
%             t=time{s1(i),s2(i),mm(i)};
            t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
            m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
            t = t_(m_==mm(i));
            Index1=0;
            for j=1:ON
                if Machine(mm(i)).GapT(j)>t                  
                    if ON==1 || j==1
                       tmp=finish{s1(i),s2(i)-1}-0; 
                    else                        
                       tmp=finish{s1(i),s2(i)-1}-Machine(mm(i)).MFT(j-1); 
                    end
                    if Machine(mm(i)).GapT(j)-t-tmp>0
                        Index1=j; %在Index工序前插入新来的工序
                        break;
                    end
                end
            end
            if Index1~=0
                Index1=Machine(mm(i)).Op(Index1);% 在s1向量中把s1(i)插入到s1(Index1)前面,Index1<i
                tmp=s1(i);
                for j=i:-1:Index1+1
                    s1(j)=s1(j-1);
                end
                s1(Index1)=tmp;
                tmp=s2(i);
                for j=i:-1:Index1+1
                    s2(j)=s2(j-1);
                end
                s2(Index1)=tmp;
                
                tmp=mm(i);
                for j=i:-1:Index1+1
                    mm(j)=mm(j-1);
                end
                mm(Index1)=tmp;
            for j=1:ON
                if Machine(mm(Index1)).Op(j)>=Index1
                    Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)+1;%由于前面插入工序，其他索引需要整体向后退一格
                end
            end
            for k=1:TM
                if k~=mm(Index1)
                    ON2=length(Machine(k).Op);
                    for h=1:ON2
                        if Machine(k).Op(h)>Index1&&Machine(k).Op(h)<i
                            Machine(k).Op(h)=Machine(k).Op(h)+1;
                        end
                    end
                end
            end
            Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
            tmp2=Machine(mm(Index1)).Op;
            Machine(mm(Index1)).Op=sort(tmp2,'ascend');
            IIndex=find(Machine(mm(Index1)).Op==Index1);
            if IIndex==1
                start{s1(Index1),s2(Index1)}=max(0,finish{s1(Index1),s2(Index1)-1});
            else
                LastOp=Machine(mm(Index1)).Op(IIndex-1);
                start{s1(Index1),s2(Index1)}=max(finish{s1(Index1),s2(Index1)-1},finish{s1(LastOp),s2(LastOp)});
            end            
            finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t;
            ON=ON+1;
            for j=1:ON
                Index1=Machine(mm(Index1)).Op(j);
                if j==1
                    Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)};                   
                else
                    LastOp=Machine(mm(Index1)).Op(j-1);
                    Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                end
                Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                    if Machine(mm(Index1)).GapT(j)<0
                        Machine(mm(Index1)).GapT(j)
                    end
            end
            mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else%如果index1==0说明没有空位就要老老实实的去在后面
                start{s1(i),s2(i)}=max(Machine(mm(i)).MFT(ON),finish{s1(i),s2(i)-1});
%                 mt{mm(i)}=start{s1(i),s2(i)}+time{s1(i),s2(i),mm(i)};   
                t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
                m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
                mt{mm(i)} = start{s1(i),s2(i)}+t_(m_==mm(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                gap=start{s1(i),s2(i)}-Machine(mm(i)).MFT(ON);
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,gap];
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];              
            end
         else
%              mt{mm(i)}=finish{s1(i),s2(i)-1}+time{s1(i),s2(i),mm(i)};
             t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
             m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
             mt{mm(i)} = finish{s1(i),s2(i)-1}+t_(m_==mm(i));
             start{s1(i),s2(i)}=finish{s1(i),s2(i)-1};
             finish{s1(i),s2(i)}=mt{mm(i)};
             Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
             Machine(mm(i)).GapT=[Machine(mm(i)).GapT,start{s1(i),s2(i)}];
             Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
         end
                       
        end
    end

    newp=s1;
end

function [newp]=AS2FAS(p_chrom,m_chrom,FJ,bkg)
% global N H TM time;
    N = bkg.job;
    H = bkg.operations';
    TM = bkg.operation;
    e=0;
    JOBN=length(FJ);
    SH=length(p_chrom);
    finish={};%工序完工时间
    start={};%工序开始时间
    for i=1:JOBN%初始化完成时间矩阵
        JOBI=FJ(i);
        for j=1:H(JOBI)
            finish{JOBI,j}=e;
            start{JOBI,j}=e;
        end
    end
    
    mt=cell(1,TM);
    for i=1:TM%初始化机器最大完成时间数组
        Machine(i).Op=[];
        Machine(i).GapT=[];
        Machine(i).MFT=[];
        mt{i}=e;
    end
    s1=p_chrom;
    s2=zeros(1,SH);
    p=zeros(1,N);
    
    for i=1:SH
        p(s1(i))=p(s1(i))+1;%记录过程是否加工完成 完成一次加一
        s2(i)=p(s1(i));%记录加工过程中，工件的次数
    end
    
    for i=1:SH
        t1=s1(i);%记录到当前是那个工件
        t2=s2(i);%记录当前工件是加工到第几次
        mm(i)=m_chrom(1,sum(H(1,1:t1-1))+t2);%提取该工序该次加工的机器选择，因为机器码的排列表示该工件第几次加工所选的机器，是一段表示一个工件
    end
    %开始解码
    for i=SH:-1:1 %倒序解码染色体
        if(s2(i)==H(s1(i))) 
         ON=length(Machine(mm(i)).Op);%该机器上目前的工序数   
         if ON>0
%             t=time{s1(i),s2(i),mm(i)};
            t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
            m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
            t = t_(m_==mm(i));
            Index1=0;
            for j=1:ON
                if Machine(mm(i)).GapT(j)-t>0
                    Index1=j; %在Index工序前插入新来的工序
                    break;
                end
            end
            if Index1~=0
                Index1=Machine(mm(i)).Op(Index1);% 在s1向量中把s1(i)插入到s1(Index1)后面,i<Index1
                tmp=s1(i);
                for j=i:Index1-1
                    s1(j)=s1(j+1);
                end
                s1(Index1)=tmp;
                tmp=s2(i);
                for j=i:Index1-1
                    s2(j)=s2(j+1);
                end
                s2(Index1)=tmp;
                
                tmp=mm(i);
                for j=i:Index1-1
                    mm(j)=mm(j+1);
                end
                mm(Index1)=tmp;
                
            for j=1:ON
                if Machine(mm(Index1)).Op(j)<=Index1
                    Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)-1;%由于前面插入工序，其他索引需要整体向后退一格
                end
            end
            for k=1:TM
                if k~=mm(Index1)
                    ON2=length(Machine(k).Op);
                    for h=1:ON2
                        if Machine(k).Op(h)<Index1&&Machine(k).Op(h)>i
                            Machine(k).Op(h)=Machine(k).Op(h)-1;
                        end
                    end
                end
            end
            Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
            tmp2=Machine(mm(Index1)).Op;
            Machine(mm(Index1)).Op=sort(tmp2,'descend');
            IIndex=find(Machine(mm(Index1)).Op==Index1);
            if IIndex==1
                start{s1(Index1),s2(Index1)}=0;
            else
                LastOp=Machine(mm(Index1)).Op(IIndex-1);
                start{s1(Index1),s2(Index1)}=max(0,finish{s1(LastOp),s2(LastOp)});
            end   
                finish{s1(Index1),s2(Index1)}=t+start{s1(Index1),s2(Index1)};
                ON=ON+1;
            for j=1:ON
                Index1=Machine(mm(Index1)).Op(j);
                if j==1
                    Machine(mm(Index1)).GapT(j)=0;                   
                else
                    LastOp=Machine(mm(Index1)).Op(j-1);
                    Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                    if Machine(mm(Index1)).GapT(j)<0
                        Machine(mm(Index1)).GapT(j)
                    end
                end
                Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
            end
            mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else %如果index1==0说明没有空位就要老老实实的去在后面
                start{s1(i),s2(i)}=Machine(mm(i)).MFT(ON);
%                 mt{mm(i)}=start{s1(i),s2(i)}+time{s1(i),s2(i),mm(i)};     
                t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
                m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
                mt{mm(i)} = start{s1(i),s2(i)}+t_(m_==mm(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
            end
            
         else
%              mt{mm(i)}=time{s1(i),s2(i),mm(i)};
             t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
             m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
             mt{mm(i)} = t_(m_==mm(i));
             start{s1(i),s2(i)}=0;
             finish{s1(i),s2(i)}=mt{mm(i)};
             Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
             Machine(mm(i)).GapT=[Machine(mm(i)).GapT,0];
             Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
         end
         
        else
            %开始在相同的机器上寻找合适的空位插入，如果自身的加工时间小于空闲时间则可以插入
         ON=length(Machine(mm(i)).Op);%该机器上目前的工序数   
         if ON>0
%             t=time{s1(i),s2(i),mm(i)};
            t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
            m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
            t = t_(m_==mm(i));
            Index1=0;
            for j=1:ON
                if Machine(mm(i)).GapT(j)>t                  
                    if ON==1 || j==1
                       tmp=finish{s1(i),s2(i)+1}-0; 
                    else                        
                       tmp=finish{s1(i),s2(i)+1}-Machine(mm(i)).MFT(j-1);
                       % tmp < 0？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
                    end
                    if Machine(mm(i)).GapT(j)-t-tmp>0
                        Index1=j; %在Index工序前插入新来的工序
                        break;
                    end
                end
            end
            if Index1~=0
                Index1=Machine(mm(i)).Op(Index1);% 在s1向量中把s1(i)插入到s1(Index1)前面,Index1<i
                tmp=s1(i);
                for j=i:Index1-1
                    s1(j)=s1(j+1);
                end
                s1(Index1)=tmp;
                tmp=s2(i);
                for j=i:Index1-1
                    s2(j)=s2(j+1);
                end
                s2(Index1)=tmp;
                
                tmp=mm(i);
                for j=i:Index1-1
                    mm(j)=mm(j+1);
                end
                mm(Index1)=tmp;
                
            for j=1:ON
                if Machine(mm(Index1)).Op(j)<=Index1
                    Machine(mm(Index1)).Op(j)=Machine(mm(Index1)).Op(j)-1;%由于前面插入工序，其他索引需要整体向后退一格
                end
            end
            for k=1:TM
                if k~=mm(Index1)
                    ON2=length(Machine(k).Op);
                    for h=1:ON2
                        if Machine(k).Op(h)<Index1&&Machine(k).Op(h)>i
                            Machine(k).Op(h)=Machine(k).Op(h)-1;
                        end
                    end
                end
            end
            Machine(mm(Index1)).Op=[Machine(mm(Index1)).Op,Index1];
            tmp2=Machine(mm(Index1)).Op;
            Machine(mm(Index1)).Op=sort(tmp2,'descend');
            IIndex=find(Machine(mm(Index1)).Op==Index1);
            if IIndex==1
                start{s1(Index1),s2(Index1)}=max(0,finish{s1(Index1),s2(Index1)+1});
            else
                LastOp=Machine(mm(Index1)).Op(IIndex-1);
                start{s1(Index1),s2(Index1)}=max(finish{s1(Index1),s2(Index1)+1},finish{s1(LastOp),s2(LastOp)});
            end            
            finish{s1(Index1),s2(Index1)}=start{s1(Index1),s2(Index1)}+t;
            ON=ON+1;
            for j=1:ON
                Index1=Machine(mm(Index1)).Op(j);
                if j==1
                    Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)};                   
                else
                    LastOp=Machine(mm(Index1)).Op(j-1);
                    Machine(mm(Index1)).GapT(j)=start{s1(Index1),s2(Index1)}-finish{s1(LastOp),s2(LastOp)};
                end
                Machine(mm(Index1)).MFT(j)=finish{s1(Index1),s2(Index1)};
                if Machine(mm(Index1)).GapT(j)<0
                    Machine(mm(Index1)).GapT(j)
                end
            end
            mt{mm(Index1)}=Machine(mm(Index1)).MFT(ON);
            else%如果index1==0说明没有空位就要老老实实的去在后面
                start{s1(i),s2(i)}=max(Machine(mm(i)).MFT(ON),finish{s1(i),s2(i)+1});
%                 mt{mm(i)}=start{s1(i),s2(i)}+time{s1(i),s2(i),mm(i)};     
                t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
                m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
                mt{mm(i)} = start{s1(i),s2(i)}+t_(m_==mm(i));
                finish{s1(i),s2(i)}=mt{mm(i)};
                Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
                gap=start{s1(i),s2(i)}-Machine(mm(i)).MFT(ON);
                Machine(mm(i)).GapT=[Machine(mm(i)).GapT,gap];
                Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];              
            end
         else
%              mt{mm(i)}=finish{s1(i),s2(i)+1}+time{s1(i),s2(i),mm(i)};
             t_ = cell2mat(bkg.times{s1(i)}(s2(i)));
             m_ = cell2mat(bkg.machines{s1(i)}(s2(i)));
             mt{mm(i)} = finish{s1(i),s2(i)+1}+t_(m_==mm(i));
             start{s1(i),s2(i)}=finish{s1(i),s2(i)+1};
             finish{s1(i),s2(i)}=mt{mm(i)};
             Machine(mm(i)).Op=[Machine(mm(i)).Op,i];
             Machine(mm(i)).GapT=[Machine(mm(i)).GapT,start{s1(i),s2(i)}];
             Machine(mm(i)).MFT=[Machine(mm(i)).MFT,mt{mm(i)}];
         end
                       
        end
    end

    newp=s1;
end
