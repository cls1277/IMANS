function [newp,newm,newf]=LS3(p_chrom,m_chrom,f_chrom,bkg)%swap Two-point exchange neigborhood
%�ڹؼ�����swap
SH = bkg.operation;
IndexO1=ceil(rand*SH);IndexO2=ceil(rand*SH);
while IndexO1==IndexO2
    IndexO2=ceil(rand*SH);
end
newp=p_chrom;
newm=m_chrom;
newf=f_chrom;
tmp=newp(IndexO1);
newp(IndexO1)=newp(IndexO2);%���������Ⱦɫ��ֵ
newp(IndexO2)=tmp;
end