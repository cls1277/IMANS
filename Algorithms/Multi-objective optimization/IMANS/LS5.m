function [newp,newm,newf]=LS5(p_chrom,m_chrom,f_chrom,fitness,bkg)
% N6
%N6�ҵ��ؼ�·���Լ�����ͬһ�������������������Ϊһ���ڽӿ飬����ͷ���β����м�飬�м���ͷβ����֮��Ĺ������ѡһ�����뵽ͷ������ǰ�����ѡһ�����뵽β�������ͷ���β��ģ���ͷ���β��������ǰ�Ĺ�����뵽β������󣬽�β��ͷ������ĺ���Ĺ�����뵽ͷ������ǰ��
N = bkg.job;
SH = bkg.operation;
    s1=p_chrom;
    s2=zeros(1,SH);
    p=zeros(1,N);
    
    newf=f_chrom;
    
    for i=1:SH
        p(s1(i))=p(s1(i))+1;%��¼�����Ƿ�ӹ���� ���һ�μ�һ
        s2(i)=p(s1(i));%��¼�ӹ������У������Ĵ���
    end
    
P1=[];P2=[];IP1=[];IP2=[];
for i=1:SH
    t1=s1(i);%��¼����ǰ���Ǹ�����
    t2=s2(i);%��¼��ǰ�����Ǽӹ����ڼ���
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

if fitness(1,3)==1
    [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(P1,m_chrom,FJ1,bkg); %�ؼ�·�����ص�������Ⱦɫ���йؼ�������±�
else
    [CriticalPath,CriticalBlock,block]=FindCriticalPathDFJSP(P2,m_chrom,FJ2,bkg);
end

for i=1:block
    BL=length(CriticalBlock(i).B);
    if BL>1
       
        if i==1 %���ѡһ��������뵽β�������
            Index1=ceil(rand*(BL-1));
            Index2=BL;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            if fitness(1,3)==1 %update P1
                tmp=P1(Index1);
                for j=Index1:Index2-1
                    P1(j)=P1(j+1);
                end
                P1(Index2)=tmp;
            else %update P2
                tmp=P2(Index1);
                for j=Index1:Index2-1
                    P2(j)=P2(j+1);
                end
                P2(Index2)=tmp;
            end
        end
        
        if i==block %���ѡһ��������뵽ͷ������ǰ
            Index1=1;
            Index2=ceil(rand*(BL-1))+1;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            if fitness(1,3)==1 %update P1
                tmp=P1(Index2);
                for j=Index2:-1:Index1+1
                    P1(j)=P1(j-1);
                end
                P1(Index1)=tmp;
            else %update P2
                tmp=P2(Index2);
                for j=Index2:-1:Index1+1
                    P2(j)=P2(j-1);
                end
                P2(Index1)=tmp;
            end
        
        end
        
        if i>1&&i<block&&BL>2 %���ѡһ��������뵽ͷ������ǰ�����ѡһ�����뵽β�������
            Index1=ceil(rand*(BL-2))+1; %�м����뵽β��
            Index2=BL;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            if fitness(1,3)==1 %update P1
                tmp=P1(Index1);
                for j=Index1:Index2-1 %β������ǰ��Ĺ�������ƣ�����ٰ�ѡ�й������β��
                    P1(j)=P1(j+1);
                end
                P1(Index2)=tmp;
            else %update P2
                tmp=P2(Index1);
                for j=Index1:Index2-1
                    P2(j)=P2(j+1);
                end
                P2(Index2)=tmp;
            end
            
            Index1=1; %�м�����ͷ��֮ǰ
            Index2=ceil(rand*(BL-2))+1;
            Index1=CriticalBlock(i).B(Index1);
            Index2=CriticalBlock(i).B(Index2);
            if fitness(1,3)==1 %update P1
                tmp=P1(Index2);
                for j=Index2:-1:Index1+1 %ͷ���������Ĺ�����ǰ���ƣ�����ٰ�ѡ�й������ͷ��
                    P1(j)=P1(j-1);
                end
                P1(Index1)=tmp;
            else %update P2
                tmp=P2(Index2);
                for j=Index2:-1:Index1+1
                    P2(j)=P2(j-1);
                end
                P2(Index1)=tmp;
            end
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