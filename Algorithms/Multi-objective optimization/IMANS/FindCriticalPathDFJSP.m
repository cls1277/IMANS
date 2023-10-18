function [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(p_chrom,m_chrom,FJ,bkg)
%�������빹��������ͼ�Ľṹ������ͼ�а����ڵ㼯 �� ���й���
%�����߼����������ڵ�֮���Ƿ��������
%ÿ���ڵ����һ���ӹ�ʱ�䡣ÿ�������ĵ�һ��������Ϊ��ʼ�ڵ����һ���ڵ㡣���ֻ��N��·����NΪ����������ÿ����������һ��������Ϊ�����ڵ����һ������
%�������Գ�����ȵ������ԣ���������ͼ��ÿ���ڵ�ֻ����������������ȡ�����ǰ��������һ�����򣬺͵�ǰ��������һ������
%���������ͼ��������ӿ�ʼ�ڵ㵽�����ڵ������깤ʱ��Ĺؼ�·����
%���õ����ݽṹΪ������������������洢ÿ��������±ꡣÿ������Ԫ�غ��һ��Ԫ��Ϊ��ǰ�����������һ�����򣬵ڶ���Ԫ��Ϊ��ǰ�����������һ������
%������Ԫ�ض�����˫���������������ݵ�Ԫ�غ��档ʵ��ΪSH*4�ľ���
%����ÿ�������ĵ�һ�����򣬴ӵ�һ������ʼ����ȱ����ڶ������һ��·��������ȱ������������һ��·����ÿ��·���Ľ�����׼���ǣ���������һ�е�ֵΪ0.

N = bkg.job;
H = bkg.operations';

JOBN=length(FJ);
SH=length(p_chrom);
% drawFJSP(p_chrom,m_chrom,FJ);

%ÿһ�зֲ���ʾ ��ǰ����������һ���������������ǰ�������ڻ�������һ��������������ǰ���򹤼�����һ���������������ǰ�������ڻ�������һ�����������
digraph=zeros(SH,4);%���ھ���������кű���Ϳ��Ա�ʾ����ֵ������ʡ��һ��
dflag=zeros(SH,4);%���ڱ�Ǹõ��Ƿ��ҵ�
    %�� ��ɽ�����ÿ������Ĺ����͹����
    s1=p_chrom;
    s2=zeros(1,SH);
    p=zeros(1,N);
    CriticalPath=[];%�ؼ�·�������а����ؼ�������±�
    
    for i=1:SH
        p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
        s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
    end
    
    for i=1:SH
        t1=s1(i);%��¼����ǰ���Ǹ�����
        t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
        mm(i)=m_chrom(1,sum(H(1,1:t1-1))+t2);%��ȡ�ù���ôμӹ��Ļ���ѡ����Ϊ����������б�ʾ�ù����ڼ��μӹ���ѡ�Ļ�������һ�α�ʾһ������
    end
    %��������ͼ
    for i=1:SH
        to=s1(i);tm=mm(i);
        %�ҵ�ǰ�����������һ������
        
        for j=i+1:SH
            if to==s1(j)&&dflag(i,1)==0
                digraph(i,1)=j;
                dflag(i,1)=1; %��ָ���Ѿ��ҵ�
                digraph(j,3)=i;
                dflag(j,3)=1;%��ָ���Ѿ��ҵ�
            end
            %�ҵ�ǰ�����������һ������
           if tm==mm(j)&&dflag(i,2)==0
                digraph(i,2)=j;
                dflag(i,2)=1;
                digraph(j,4)=i;
                dflag(j,4)=1;
           end
            %ȫ��������������ǰѭ�� ���������ǰ����ǰ�����������һ������
           if (dflag(i,1)==1||s2(i)==H(s1(i)))&&dflag(i,2)==1
               break;
           end
           
        end
    end
    VELTime=cell(SH,2);%ÿ����������翪ʼʱ�������ʼʱ��
    e=0;
    %���ù�ȱ�������ÿ����������翪ʼʱ�������ʼʱ��,��������ù�����ȱ���������ڸ��º������ʱ�򣬳���ǰ�������������㣬�����¼������
    %������ȱ�����˳���ǣ��������й����ĵ�һ�������������й����ĵڶ�������
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
    
    %�������翪ʼʱ��
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
       if t1>t2 %���t1����t2��
           VELTime{Index,1}=t1;
       else
           VELTime{Index,1}=t2;
       end
       
       if digraph(Index,1)==0 %����ǹ��������һ��������˵����һ���ڵ���LastNode��������ڵ�, ������ڵ����һ���ڵ������й��������һ������
%            t1=time{s1(Index),s2(Index),mm(Index)}+VELTime{Index,1};%ѡ������Ϊ��������ڵ�
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            t1 = t_(m_==mm(Index)) + VELTime{Index,1};
           if t1>LastNode
               LastNode=t1;
           end
       end
    end
    %��������ʼʱ��
    for i=SH:-1:1
        Index=level(i); 
%         fprintf('%s %d %s %d\r\n','O',s1(Index),'.',s2(Index));
        if digraph(Index,1)==0
            if digraph(Index,2)==0
%                 VELTime{Index,2}=LastNode-time{s1(Index),s2(Index),mm(Index)}; %���һ����������Ľ��Ĺ��������ʼʱ��������翪ʼʱ��             
            t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
            m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
            VELTime{Index,2}=LastNode-t_(m_==mm(Index)); %���һ����������Ľ��Ĺ��������ʼʱ��������翪ʼʱ��
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
                if t1>t2 %ѡ��С��
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
%            t1=VELTime{next,2}-time{s1(Index),s2(Index),mm(Index)};%��һ���ڵ����ǰ�����ߵ�ʱ�䣬����ʱ���ڽڵ������м���ǰ�ڵ���ʱ��
          t_ = cell2mat(bkg.times{s1(Index)}(s2(Index)));
          m_ = cell2mat(bkg.machines{s1(Index)}(s2(Index)));
          t1=VELTime{next,2}-t_(m_==mm(Index));%��һ���ڵ����ǰ�����ߵ�ʱ�䣬����ʱ���ڽڵ������м���ǰ�ڵ���ʱ��
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
       if t1>t2 %ѡ��С��
           VELTime{Index,2}=t2;
       else
           VELTime{Index,2}=t1;
       end
%        if VELTime{Index,2}<0
%            fprintf('%s %d %s %d\r\n','O',s1(Index),'.',s2(Index));
%        end
    end
    
    %������ʼʱ������翪ʼʱ�䣬ʣ��ʱ����С�Ľ�㣬��Ϊ�ؼ�·���ϵĽ�㡣
    Idletime=cell(SH,1);
    for i=1:SH
        Idletime{i,1}=VELTime{i,2}-VELTime{i,1};
    end
    %�ҵ��ɳ�ʱ����С�Ĺ�����Ϊ�ؼ�·���ϵĹ���
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
    CriticalBlock(block).B=[];%��ʾ��ͬһ�������Ͼ���ǰ���ϵ��˳��� ����
    CriticalBlock(block).B=[CriticalBlock(block).B,CriticalPath(1)];m0=mm(CriticalPath(1));
    for i=2:L
        m=mm(CriticalPath(i));
        if m0==m
            CriticalBlock(block).B=[CriticalBlock(block).B,CriticalPath(i)]; %��ͬһ�����������ڵĹ������Ϊһ���ڽӿ�
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
       if t1>t2 %���t1����t2��
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
       if t1>t2 %���t1����t2��
           VELTime{Index,1}=t2;
       else
           VELTime{Index,1}=t1;
       end
       value=VELTime{Index,1};
    
end
