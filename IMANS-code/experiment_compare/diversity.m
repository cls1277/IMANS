function DV=diversity(obj, fmin, fmax)%Ϊ�˶Աȷ���ҲҪ��һ��
    [N,~]  = size(obj);
    obj = (obj-repmat(fmin,N,1))./repmat(fmax-fmin,N,1);
    obj_size = N;
    distance=zeros(obj_size-1,1);%ֻ�����ڵľ���
    d_averge=0;
    DV=0;
    if obj_size==1
        DV = 1;
        return ;
    end
    A=obj(1,:); B=[1 0];
    df = sqrt((A-B)*(A-B)');
    A=obj(obj_size,:); B=[0 1];
    dl = sqrt((A-B)*(A-B)');
%     df = 0; dl = 0;
for i=1:obj_size-1
     A=obj(i,:);B=obj(i+1,:);
     distance(i,1)=(A-B)*(A-B)';             % ����ŷʽ���룬����û�п����ţ���Ϊֻ����Ҫ����С����һ����Ч��
     distance(i,1)=sqrt(distance(i,1));
     d_averge=d_averge+distance(i,1);
end
d_averge=d_averge/(obj_size-1);
for i=1:obj_size-1
    DV=DV+abs(distance(i,1)-d_averge);
end
% DV=DV/((obj_size-1)*d_averge);
DV = (DV+df+dl)/((obj_size-1)*d_averge+df+dl);
if obj_size<1 %��һ���ǳ�С����ȡ��0
    DV=0.00001;
end
end