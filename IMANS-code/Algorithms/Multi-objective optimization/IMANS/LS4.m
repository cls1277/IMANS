function [newp,newm,newf]=LS4(p_chrom,m_chrom,f_chrom, bkg)%one point insert, random insert neighborhood
%�ڷǹؼ�����insert
SH = bkg.operation;

IndexO1=ceil(rand*SH);IndexO2=ceil(rand*SH);
while IndexO1==IndexO2
    IndexO2=ceil(rand*SH);
end

newp=p_chrom;
newm=m_chrom;
newf=f_chrom;

if IndexO1>IndexO2
    tmp=IndexO1;
    IndexO1=IndexO2;
    IndexO2=tmp;
end

tmp=newp(IndexO2);%ȡ��t2���ֵ�ֵ
%��t1+1��t2����ֵ�������
for i=IndexO2:-1:IndexO1+1
    newp(i)=newp(i-1);
end
newp(IndexO1)=tmp;%��t2��ֵ���뵽t1��

end