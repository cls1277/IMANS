function [newp,newm,newf]=LS3(p_chrom,m_chrom,f_chrom,bkg)%swap Two-point exchange neigborhood
%在关键工厂swap
SH = bkg.operation;
IndexO1=ceil(rand*SH);IndexO2=ceil(rand*SH);
while IndexO1==IndexO2
    IndexO2=ceil(rand*SH);
end
newp=p_chrom;
newm=m_chrom;
newf=f_chrom;
tmp=newp(IndexO1);
newp(IndexO1)=newp(IndexO2);%交换两点的染色体值
newp(IndexO2)=tmp;
end